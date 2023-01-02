from __future__ import annotations

import os
import textwrap
import warnings
from functools import partial
from numbers import Number
from typing import Optional, Any, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import torch
from hdxrate import k_int_from_sequence
from scipy import constants
from scipy.constants import R
from scipy.integrate import solve_ivp

from pyhdx.alignment import align_dataframes
from pyhdx.fileIO import dataframe_to_file
from pyhdx.process import verify_sequence, parse_temperature
from pyhdx.support import reduce_inter, dataframe_intersection, array_intersection
from pyhdx.config import cfg


def protein_wrapper(func, *args, **kwargs):
    warnings.warn("Protein wrapper and objects are deprecated", DeprecationWarning)
    metadata = kwargs.pop("metadata", {})
    [metadata.update(arg.metadata) for arg in args if isinstance(arg, Protein)]

    df_args = [arg.df if isinstance(arg, Protein) else arg for arg in args]
    return_value = func(*df_args, **kwargs)
    if isinstance(return_value, pd.DataFrame):
        return Protein(return_value, **metadata)

    return return_value


class Protein(object):
    """Object describing a protein

    Protein objects are based on panda's DataFrame's with added functionality

    Parameters
    ----------
    data : :class:`~numpy.ndarray` or :obj:`dict` or :class:`~pandas.DataFrame`
        data object to initiate the protein object from
    index : :obj:`str`, optional
        Name of the column with the residue number (index column)

    **metadata
        Dictionary of optional metadata.


    """

    def __init__(self, data, index=None, **metadata):
        warnings.warn("Protein wrapper and objects are deprecated", DeprecationWarning)
        self.metadata = metadata
        if isinstance(data, dict) or isinstance(data, np.ndarray):
            self.df = pd.DataFrame(data)
            self.df.set_index(index, inplace=True)
        elif isinstance(data, pd.DataFrame):
            self.df = data.copy()
            if not self.df.index.is_integer():
                raise ValueError(
                    f"Invalid index type {type(self.df.index)} for supplied DataFrame, must be integer index"
                )

        if not self.df.index.is_unique:
            raise ValueError("Protein dataframe indices must be unique")

        new_index = pd.RangeIndex(
            start=self.df.index.min(), stop=self.df.index.max() + 1, name="r_number"
        )
        self.df = self.df.reindex(new_index)

    def __str__(self):
        s = self.df.__str__()
        try:
            full_s = f"Protein {self.metadata['name']}\n" + s
            return full_s
        except KeyError:
            return s

    def __len__(self):
        return len(self.df)

    def __getattr__(self, item):
        attr = getattr(self.df, item)
        if callable(attr):
            return partial(protein_wrapper, attr, metadata=self.metadata)
        else:
            return attr

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__.update(d)

    def _make_protein(self, df_out, other):
        """Make a new :class:`~pyhdx.models.Protein` object and combine metadata with other metadata"""
        metadata = {**self.metadata, **other.metadata}
        protein_out = Protein(df_out, index=df_out.index.name, **metadata)
        return protein_out

    def to_file(
        self,
        file_path,
        include_version=True,
        include_metadata=True,
        fmt="csv",
        **kwargs,
    ):
        """
        Write Protein data to file.


        Parameters
        ----------
        file_path : :obj:`str`
            File path to create and write to.
        include_version : :obj:`bool`
            Set ``True`` to include PyHDX version and current time/date
        fmt : :obj:`str`
            Formatting to use, options are 'csv' or 'pprint'
        include_metadata : :obj:`bool`
            If `True`, the objects' metadata is included
        **kwargs : :obj:`dict`, optional
            Optional additional keyword arguments passed to `df.to_csv`
        Returns
        -------
        None

        """

        metadata = self.metadata if include_metadata else include_metadata
        dataframe_to_file(
            file_path,
            self.df,
            include_version=include_version,
            include_metadata=metadata,
            fmt=fmt,
            **kwargs,
        )

    def get_k_int(self, temperature, pH, **kwargs):
        """
        Calculates the intrinsic rate of the sequence. Values of no coverage or prolines are assigned a value of -1
        The rates run are for the first residue (1) up to the last residue that is covered by peptides

        When the previous residue is unknown the current residue is also assigned a value of -1.g

        Parameters
        ----------
        temperature : :obj:`float`
            Temperature of the labelling reaction (Kelvin)
        pH : :obj:`float`
            pH of the labelling reaction

        Returns
        -------

        k_int : :class:`~numpy.ndarray`
            Array of intrisic exchange rates

        """

        if "sequence" not in self:
            raise ValueError(
                "No sequence data available to calculate intrinsic exchange rates."
            )

        sequence = list(
            self["sequence"]
        )  # Includes 'X' padding at cterm if cterm > last peptide
        k_int_array = k_int_from_sequence(sequence, temperature, pH, **kwargs)
        k_int = pd.Series(k_int_array, index=self.df.index)

        return k_int

    @property
    def c_term(self):
        return self.df.index.max()

    @property
    def n_term(self):
        return self.df.index.min()

    def __getitem__(self, item):
        return self.df.__getitem__(item)

    def __setitem__(self, index, value):
        self.df.__setitem__(index, value)

    def __contains__(self, item):
        return self.df.__contains__(item)

    def __sub__(self, other):
        return protein_wrapper(self.df.subtract, other, metadata=self.metadata)

    def __add__(self, other):
        return protein_wrapper(self.df.add, other, metadata=self.metadata)

    def __truediv__(self, other):
        return protein_wrapper(self.df.truediv, other, metadata=self.metadata)

    def __floordiv__(self, other):
        return protein_wrapper(self.df.floordiv, other, metadata=self.metadata)

    def __mul__(self, other):
        return protein_wrapper(self.df.mul, other, metadata=self.metadata)


class Coverage(object):
    """
    Object describing layout and coverage of peptides and generating the corresponding matrices.
    Peptides should all belong to the same state and have the same exposure time.

    Args:
        data: DataFrame with input peptides
        weight_exponent: Value of the exponent for weighted averaging used in converting
            peptide RFU values to residue-level RFU values. Default value is 1., corresponding
            to weighted averaging where weights are the inverse of peptide length. Generally,
            weights are 1/(peptide_length)**exponent.
        n_term: Residue index of the N-terminal residue. Default value is 1, can be
            negative to accommodate for N-terminal purification tags
        c_term: Residue index number of the C-terminal residue (where first residue has
            index number 1)
        sequence: Amino acid sequence of the protein in one-letter FASTA encoding.
            Optional, if not specified the amino acid sequence from the peptide data is used
            to (partially) reconstruct the sequence. Supplied amino acid sequence must be
            compatible with sequence information in the peptides.

    """

    X: np.ndarray
    """
    Np x Nr matrix (peptides x residues). Values are 1 where residue j is in peptide i.
    """

    Z: np.ndarray
    """
    Np x Nr matrix (peptides x residues). Values are 1/(ex_residues) where residue j 
    is in peptide i.
    """
    # todo account for prolines: so that rows sum to 1 is currently not true

    def __init__(
        self,
        data: pd.DataFrame,
        n_term: Optional[int] = None,
        c_term: Optional[int] = None,
        sequence: Optional[str] = None,
    ) -> None:

        for field in ["exposure", "state"]:
            if field in data and len(np.unique(data["exposure"])) != 1:
                raise ValueError(f"Entries in field {field!r} must be unique")
        try:
            self.data = data.sort_values(["_start", "_stop"], axis=0)
        except KeyError:
            self.data = data.sort_values(["start", "stop"], axis=0)
        self.data.index.name = "peptide_id"  # todo check these are the same as parent object peptide_id (todo make wide instead of instersection)

        seq_full, seq_r = verify_sequence(data, sequence, n_term, c_term)

        # todo check if this is always correctly determined (n terminal residues usw)
        exchanges = [
            s.isupper() and (s != "X") for s in seq_r
        ]  # Boolean array True if residue exchanges, full length
        coverage = seq_r != "X"  # Boolean array for coverage
        protein_df = pd.DataFrame(
            {"sequence": seq_full, "coverage": coverage, "exchanges": exchanges},
            index=seq_full.index,
        )

        # Inclusive, exclusive interval of peptides coverage across the whole protein
        self.interval = (np.min(self.data["_start"]), np.max(self.data["_stop"]))
        self.protein = Protein(protein_df, index="r_number")

        # matrix dimensions N_peptides N_residues, dtype for PyTorch compatibility
        _exchanges = self["exchanges"]  # Array only on covered part
        self.X = np.zeros(
            (len(self.data), self.interval[1] - self.interval[0]), dtype=int
        )
        self.Z = np.zeros_like(self.X, dtype=float)
        for row, idx in enumerate(self.data.index):
            # start, end are already corrected for drop_first parameter
            start, end = self.data.loc[idx, "_start"], self.data.loc[idx, "_stop"]
            i0, i1 = self.r_number.get_loc(start), self.r_number.get_loc(end - 1)
            # i0, i1 = np.searchsorted(self.r_number, (entry['start'], entry['end']))
            self.X[row][i0 : i1 + 1] = 1
            self.Z[row][i0 : i1 + 1] = _exchanges[i0 : i1 + 1]
        self.Z = self.Z / self.data["ex_residues"].to_numpy()[:, np.newaxis]

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, item) -> pd.Series:
        """Gets columns from underlying protein and crops to interval.

        Crop interval is equal to the coverage range of peptides in this :class:`.Coverage`
        object.

        """
        pd_series = self.protein[item]
        return self.apply_interval(pd_series)

    def apply_interval(
        self, array_or_series: Union[np.ndarray, pd.Series]
    ) -> pd.Series:
        """Returns the section of `array_or_series` in the interval


        Given a Numpy array or Pandas series with a length equal to the full protein,
        returns the section of the array equal to the covered
        region. Returned series length is equal to number of columns in the X matrix

        Args:
            array_or_series: Input data object to crop to interval

        Returns:
            Input object cropped to interval of the interval spanned by the peptides

        """

        if isinstance(array_or_series, np.ndarray):
            series = pd.Series(array_or_series, index=self.protein.df.index)
            assert len(array_or_series) == len(self.protein)
        else:
            series = array_or_series

        # - 1 because interval is inclusive, exclusive and .loc slices inclusive, inclusive
        covered_slice = series.loc[self.interval[0] : self.interval[1] - 1]

        return covered_slice

    @property
    def percent_coverage(self) -> float:
        """Percentage of residues covered by peptides"""
        return 100 * np.mean(self.protein["coverage"])

    @property
    def redundancy(self) -> float:
        """Average redundancy of peptides in regions with at least 1 peptide"""
        x_coverage = self.X[:, self["coverage"]]
        return float(np.mean(np.sum(x_coverage, axis=0)))

    @property
    def avg_peptide_length(self) -> float:
        """Average length of the peptides"""
        return (self.data["end"] - self.data["start"]).mean()

    @property
    def Np(self) -> int:
        """Number of peptides."""
        return self.X.shape[0]

    @property
    def Nr(self) -> int:
        """Total number of residues spanned by the peptides."""
        return self.X.shape[1]

    # TODO homogenize this and next property
    @property
    def r_number(self) -> pd.RangeIndex:
        """Pandas index numbers corresponding to the part of the protein covered by peptides"""
        return pd.RangeIndex(self.interval[0], self.interval[1], name="r_number")

    @property
    def index(self) -> pd.RangeIndex:
        """Pandas index numbers corresponding to the part of the protein covered by peptides"""
        return self.r_number

    @property
    def block_length(self) -> np.ndarray:
        """:class:`~numpy.ndarary`: Lengths of unique blocks of residues in the peptides map,
        along the `r_number` axis"""

        # indices are start and stop values of blocks
        indices = np.sort(np.concatenate([self.data["_start"], self.data["_stop"]]))
        # indices of insertion into r_number vector gives us blocks with taking prolines into account.
        diffs = np.diff(np.searchsorted(self.r_number, indices))

        block_length = diffs[diffs != 0]
        return block_length

    @property
    def X_norm(self) -> np.ndarray:
        """:class:`~numpy.ndarray`: `X` coefficient matrix normalized column wise."""
        return self.X / np.sum(self.X, axis=0)[np.newaxis, :]

    @property
    def Z_norm(self) -> np.ndarray:
        """:class:`~numpy.ndarray`: `Z` coefficient matrix normalized column wise."""
        wts = self.Z**cfg.analysis.weight_exponent
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            z_norm = wts / np.sum(wts, axis=0)[np.newaxis, :]

        return z_norm

    def get_sections(self, gap_size: int = -1) -> list[tuple[int, int]]:
        """Get the intervals of independent sections of coverage.

        Intervals are inclusive, exclusive.
        Gaps are defined with `gap_size`, adjacent peptides with distances bigger than this value are considered not to
        overlap. Set to -1 to treat touching peptides as belonging to the same section.

        Args:
            gap_size: The size which defines a gap


        """
        intervals = [(s, e) for s, e in zip(self.data["_start"], self.data["_stop"])]
        sections = reduce_inter(intervals, gap_size=gap_size)

        return sections


class HDXMeasurement(object):
    """Main HDX data object.

    This object has peptide data of a single state and with multiple timepoints.
    Timepoint data is split into :class:`~pyhdx.models.PeptideMeasurements` objects for
    each timepoint. Supplied data is made 'uniform' such that all timepoints have the same peptides

    Args:
        data: Dataframe with all peptides belonging to a single state.
        **metadata: Dictionary of optional metadata. By default, holds the `temperature` and `pH` parameters.

    """

    # TODO alphabetize?
    data: pd.DataFrame
    """Dataframe with all peptides"""

    state: str
    """Protein state label for this HDX measurement"""

    timepoints: np.ndarray
    """Deuterium exposure times"""

    peptides: list[HDXTimepoint]
    """List of :class:`.HDXTimepoint`, one per exposure timepoint"""

    coverage: Coverage
    """Coverage object describing peptide layout"""

    def __init__(self, data: pd.DataFrame, **metadata: Any):
        self.metadata = metadata
        assert len(data["state"].unique()) == 1
        self.state = str(data["state"].iloc[0])
        self.timepoints = np.sort(np.unique(data["exposure"]))

        # todo sort happens twice now
        data = data.sort_values(["start", "stop", "sequence", "exposure"])

        # Obtain the intersection of peptides per timepoint
        df_list = [(data[data["exposure"] == exposure]) for exposure in self.timepoints]

        intersected_data = dataframe_intersection(df_list, by=["start", "stop"])

        cov_kwargs = {
            kwarg: metadata.get(kwarg) for kwarg in ["c_term", "n_term", "sequence"]
        }
        self.peptides = [HDXTimepoint(df, **cov_kwargs) for df in intersected_data]

        # Create coverage object from the first time point (as all are now equal)
        self.coverage = Coverage(intersected_data[0], **cov_kwargs)

        if self.temperature and self.pH:
            k_int = self.coverage.protein.get_k_int(self.temperature, self.pH)
            self.coverage.protein["k_int"] = k_int

        self.data = pd.concat(intersected_data, axis=0, ignore_index=True).sort_values(
            ["start", "stop", "sequence", "exposure"]
        )
        self.data["peptide_id"] = self.data.index % self.Np
        self.data.index.name = (
            "peptide_index"  # index is original index which continues along exposures
        )
        self.data_wide = (
            self.data.pivot(index="peptide_id", columns=["exposure"])
            .reorder_levels([1, 0], axis=1)
            .sort_index(axis=1, level=0, sort_remaining=False)
        )

    def __str__(self) -> str:
        """String representation of this HDX measurement object.

        Returns:
            Multiline string describing this HDX Measurement object

        """

        timepoints = ", ".join([f"{t:.2f}" for t in self.timepoints])

        s = f"""
        HDX Measurement:     {self.name}
        
        Number of peptides:         {self.Np}
        Number of residues:         {self.Nr} ({self.coverage.interval[0]} - {self.coverage.interval[1]})
        Number of timepoints:       {self.Nt}
        Timepoints:                 {timepoints} seconds
        Coverage Percentage:        {self.coverage.percent_coverage:.2f}
        Average redundancy:         {self.coverage.redundancy:.2f}   
        Average peptide length:     {self.coverage.avg_peptide_length:.2f}
        Repeatability (mean std):   {self.data['uptake_sd'].mean():.2f} Da
        Temperature:                {self.temperature} K
        pH:                         {self.pH}             
        """

        return textwrap.dedent(s.lstrip("\n"))

    def _repr_markdown_(self) -> str:
        """Markdown representation this HDX measurement object"""
        s = str(self)
        s = s.replace("\n", "<br>")
        return s

    @property
    def name(self) -> str:
        """HDX Measurement name"""
        return self.metadata.get("name", self.state)

    @property
    def temperature(self) -> Optional[float]:
        """Temperature of the H/D exchange reaction (K)."""
        temperature = self.metadata.get("temperature")
        if isinstance(temperature, (Number, type(None))):
            return temperature
        elif isinstance(temperature, dict):
            return parse_temperature(**temperature)

        return self.metadata.get("temperature", None)

    @property
    def pH(self) -> Optional[float]:
        """pH of the H/D exchange reaction."""
        return self.metadata.get("pH", None)

    @property
    def Np(self) -> int:
        """Number of peptides."""
        return self.coverage.Np

    @property
    def Nr(self) -> int:
        """Total number of residues spanned by the peptides."""
        return self.coverage.Nr

    @property
    def Nt(self) -> int:
        """Number of timepoints."""
        return len(self.timepoints)

    def __len__(self) -> int:
        import warnings

        warnings.warn("Use hdxm.Nt instead", DeprecationWarning)
        return len(self.timepoints)

    def __iter__(self):
        return self.peptides.__iter__()

    def __getitem__(self, item):
        return self.peptides.__getitem__(item)

    @property
    def rfu_residues(self) -> pd.DataFrame:
        """Relative fractional uptake per residue.

        Shape of the returned DataFrame is Nr (rows) x Nt (columns)
        """
        df = pd.concat([v.rfu_residues for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"

        return df

    @property
    def rfu_residues_sd(self) -> pd.DataFrame:
        """Standard deviations of relative fractional uptake per residue.

        Shape of the returned DataFrame is Nr (rows) x Nt (columns)
        """

        df = pd.concat([v.rfu_residues_sd for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"

        return df

    @property
    def rfu_peptides(self) -> pd.DataFrame:
        """Relative fractional uptake per peptide.

        Shape of the returned DataFrame is Np (rows) x Nt (columns)
        """
        df = pd.concat([v.rfu_peptides for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"
        return df

    @property
    def d_exp(self) -> pd.DataFrame:
        """D-uptake values (corrected for back-exchange).

        Shape of the returned DataFrame is Np (rows) x Nt (columns)
        """
        df = pd.concat([v.d_exp for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"
        return df

    # todo check shapes of k_int and timepoints, compared to their shapes in hdxmeasurementset
    def get_tensors(
        self, exchanges: bool = False, dtype: Optional[torch.dtype] = None
    ) -> dict[str, torch.Tensor]:
        """Returns a dictionary of tensor variables for fitting HD kinetics.

        Tensor variables are (shape):
        Temperature (1 x 1)
        X (Np x Nr)
        k_int (Nr x 1)
        timepoints (1 x Nt)
        d_exp (D) (Np x Nt)

        Args:
            exchanges: If ``True`` only returns tensor data describing residues which exchange
                (ie have peptides and are not prolines)
            dtype: Optional Torch data type. Use torch.float32 for faster fitting of large data
                sets, possibly at the expense of accuracy

        Returns:
            Dictionary with tensors

        """

        if "k_int" not in self.coverage.protein:
            raise ValueError(
                "Unknown intrinsic rates of exchange, please supply pH and temperature parameters"
            )
        try:
            d_exp = self.d_exp  # noqa
        except ValueError:
            raise ValueError("HDX data is not corrected for back exchange.")

        if exchanges:
            # this could be a method on coverage object similar to apply_interval; select exchanging
            bools = self.coverage["exchanges"].to_numpy()
        else:
            bools = np.ones(self.Nr, dtype=bool)

        dtype = dtype or cfg.TORCH_DTYPE
        device = cfg.TORCH_DEVICE

        tensors = {
            "temperature": torch.tensor(
                [self.temperature], dtype=dtype, device=device
            ).unsqueeze(-1),
            "X": torch.tensor(self.coverage.X[:, bools], dtype=dtype, device=device),
            "k_int": torch.tensor(
                self.coverage["k_int"].to_numpy()[bools], dtype=dtype, device=device
            ).unsqueeze(-1),
            "timepoints": torch.tensor(
                self.timepoints, dtype=dtype, device=device
            ).unsqueeze(0),
            "d_exp": torch.tensor(self.d_exp.to_numpy(), dtype=dtype, device=device),
        }

        return tensors

    def guess_deltaG(self, rates: pd.Series, correct_c_term: bool = True) -> pd.Series:
        """Obtain ΔG initial guesses from apparent H/D exchange rates.

        Units of  rates are per second.
        As the intrinsic rate of exchange of the c-terminal residue is ~100 fold lower,
        guess values for PF and ΔG are also much lower. Use the option `correct_c_term` to
        set the c-terminal guess value equal to the value of the residue preceding it.

        Args:
            rates: Apparent exchange rates (units s^-1). Series index is protein residue number.
            correct_c_term: If ``True``, sets the guess value of the c-terminal residue to the
                value of the residue preceding it.

        Returns:
            ΔG guess values (units kJ/mol)

        """
        if "k_int" not in self.coverage.protein:
            raise ValueError(
                "Unknown intrinsic rates of exchange, please supply pH and temperature parameters"
            )
        if not isinstance(rates, pd.Series):
            raise TypeError("Rates input type should be pandas.Series")

        p_guess = (self.coverage.protein["k_int"] / rates) - 1

        p_guess.clip(
            0.0, None, inplace=True
        )  # Some initial guesses might have negative PF values
        with np.errstate(divide="ignore"):
            deltaG = np.log(p_guess) * constants.R * self.temperature

        deltaG.replace([np.inf, -np.inf], np.nan, inplace=True)

        if correct_c_term and self.coverage.protein.c_term in deltaG.index:
            deltaG.loc[self.coverage.protein.c_term] = deltaG.loc[
                self.coverage.protein.c_term - 1
            ]

        return deltaG

    def to_file(
        self,
        file_path: os.PathLike,
        include_version: bool = True,
        include_metadata: bool = True,
        fmt: str = "csv",
        **kwargs: Any,
    ) -> None:
        """Write the data in this :class:`.HDXMeasurement` to file.

        Args:
            file_path: File path to create and write to.
            include_version: Set ``True`` to include PyHDX version and current time/date
            fmt: Formatting to use, options are 'csv' or 'pprint'
            include_metadata: If ``True``, the objects' metadata is included
            **kwargs: Optional additional keyword arguments passed to `df.to_csv`

        """

        # requires testing dont think this works as intended
        # should use self.metadata if include_metadata is the bool True otherwise if its a dict use that
        metadata = self.metadata if include_metadata else include_metadata
        df = self.data
        dataframe_to_file(
            file_path,
            df,
            include_version=include_version,
            include_metadata=metadata,
            fmt=fmt,
            **kwargs,
        )


class HDXTimepoint(Coverage):
    """Class with subset of peptides corresponding to only one state and exposure

    Args:
        data: Dataframe with input data
        **kwargs:

    """

    state: str
    """Protein state label for this HDX timepoint"""

    exposure: float
    """Deuterium exposure time for this HDX timepoint (units seconds)"""

    def __init__(self, data: pd.DataFrame, **kwargs: Any) -> None:
        assert len(np.unique(data["exposure"])) == 1, "Exposure entries are not unique"
        assert len(np.unique(data["state"])) == 1, "State entries are not unique"

        super(HDXTimepoint, self).__init__(data, **kwargs)

        self.state = self.data["state"][0]
        self.exposure = self.data["exposure"][0]

    @property
    def rfu_peptides(self) -> pd.Series:
        """Relative fractional uptake per peptide"""
        return self.data["rfu"]

    @property
    def d_exp(self) -> pd.Series:
        """Experimentally measured D-values (corrected)"""
        return self.data["uptake_corrected"]

    @property
    def name(self) -> str:
        """Name of this HDX timepoint

        Format is <state>_<exposure>
        """
        return f"{self.state}_{self.exposure}"

    @property
    def rfu_residues(self) -> pd.Series:
        """Relative fractional uptake (RFU) per residue.

        RFU values are obtained by weighted averaging, weight value is the length of
        each peptide

        """
        return self.weighted_average("rfu")

    @property
    def rfu_residues_sd(self) -> pd.Series:
        """Error propagated standard deviations of RFU per residue."""

        return self.propagate_errors("rfu sd")

    # todo allow pd.Series?
    def calc_rfu(self, residue_rfu: np.ndarray) -> np.ndarray:
        """
        Calculates RFU per peptide given an array of individual residue scores

        Parameters
        ----------
        residue_rfu : :class:`~numpy.ndarray`
            Array of rfu per residue of length `prot_len`

        Returns
        -------

        rfu : :class:`~numpy.ndarray`
            Array of rfu per peptide
        """

        rfu = self.Z.dot(residue_rfu)
        return rfu

    def weighted_average(self, field: str) -> pd.Series:
        """
        Calculate per-residue weighted average of values in data column

        Args:
            field: Data field (column) to calculated weighted average of

        Returns:
            THe weighted averaging result

        """

        array = self.Z_norm.T.dot(self.data[field])
        series = pd.Series(array, index=self.index)

        return series

    def propagate_errors(self, field: str) -> pd.Series:
        """Propagate errors on `field` when calculating per-residue weighted average values.

        Args:
            field: Data field (column) of errors to propagate.

        Returns:
            Propagated errors per residue.

        """

        array = np.sqrt((self.Z_norm**2).T.dot(self.data[field] ** 2))
        series = pd.Series(array, index=self.index)

        return series


class CoverageSet(object):
    """Coverage object for multiple :class:`.HDXMeasurement` objects.

    This objects finds the minimal interval of residue numbers which fit all :class:`.HDXMeasurement`s


    Args:
        hdxm_list: List of input :class:`.HDXMeasurment objects.

    """

    # todo perhaps this object should have X
    def __init__(self, hdxm_list: list[HDXMeasurement]):
        self.hdxm_list = hdxm_list

        # todo create Coverage object for the 3d case
        intervals = np.array(
            [hdxm_list.coverage.interval for hdxm_list in self.hdxm_list]
        )
        self.interval = (intervals[:, 0].min(), intervals[:, 1].max())
        # TODO should be pandas dataframe? or should be the same on both coverage objects
        self.r_number = np.arange(*self.interval)

        # TODO properties?
        self.Ns = len(self.hdxm_list)
        self.Nr = len(self.r_number)
        self.Np = np.max([hdxm.Np for hdxm in self.hdxm_list])
        self.Nt = np.max([hdxm.Nt for hdxm in self.hdxm_list])

    # TODO in subclass
    @property
    def index(self) -> pd.RangeIndex:
        """Index of residue numbers"""
        return pd.RangeIndex(self.interval[0], self.interval[1], name="r_number")

    # TODO in subclass
    def apply_interval(self, array_or_series):
        """Given a Numpy array or Pandas series with a length equal to the full protein, returns the section of the array equal to the covered
        region. Returned series length is equal to number of columns in the X matrix

        """
        # todo testing and 2d array support
        if isinstance(array_or_series, np.ndarray):
            series = pd.Series(array_or_series, index=self.index)
            assert len(array_or_series) == len(self.index)
        else:
            series = array_or_series

        # - 1 because interval is inclusive, exclusive and .loc slices inclusive, inclusive
        covered_slice = series.loc[self.interval[0] : self.interval[1] - 1]

        return covered_slice

    @property
    def s_r_mask(self) -> np.ndarray:
        """Sample-residue mask

        Boolean array where entries `ij` are ``True`` if residue `j` is covered by peptides of
        sample `i` (Coverage aps not taken into account)
        """

        mask = np.zeros((self.Ns, self.Nr), dtype=bool)
        for i, hdxm in enumerate(self.hdxm_list):
            interval_sample = hdxm.coverage.interval
            i0 = interval_sample[0] - self.interval[0]
            i1 = interval_sample[1] - self.interval[0]

            mask[i, i0:i1] = True

        return mask

    def get_masks(self) -> dict[str, np.ndarray]:
        """mask of shape NsxNr with True entries covered by hdx measurements (exluding gaps)"""
        sr_mask = np.zeros((self.Ns, self.Nr), dtype=bool)
        st_mask = np.zeros((self.Ns, self.Nt), dtype=bool)
        spr_mask = np.zeros((self.Ns, self.Np, self.Nr), dtype=bool)
        spt_mask = np.zeros((self.Ns, self.Np, self.Nt), dtype=bool)
        for i, hdxm in enumerate(self.hdxm_list):
            interval_sample = hdxm.coverage.interval
            i0 = interval_sample[0] - self.interval[0]
            i1 = interval_sample[1] - self.interval[0]

            sr_mask[i, i0:i1] = True
            st_mask[i, -hdxm.Nt :] = True
            spr_mask[i, 0 : hdxm.Np, i0:i1] = True
            spt_mask[i, 0 : hdxm.Np, -hdxm.Nt :] = True

        mask_dict = {"sr": sr_mask, "st": st_mask, "spr": spr_mask, "spt": spt_mask}

        return mask_dict


class HDXMeasurementSet(object):
    """
    Set of multiple :class:`~pyhdx.models.HDXMeasurement` s

    Args:
        hdxm_list: Input list of :class:`.HDXMeasurement`

    """

    timepoints: np.ndarray
    """
    Array with timepoints, shape is Ns x Nt, padded with zeros in case of samples with
    unequal number of timepoints
    """

    # TODO this is a property on HDXMeasurement
    d_exp: np.ndarray
    """
    Array with measured D-uptake values, shape is Ns x Np x Nt, padded with zeros.
    """

    def __init__(self, hdxm_list: list[HDXMeasurement]) -> None:
        self.hdxm_list = hdxm_list

        self.coverage = CoverageSet(hdxm_list)
        self.masks = self.coverage.get_masks()

        timepoints_values = np.concatenate([hdxm.timepoints for hdxm in self.hdxm_list])
        self.timepoints = np.zeros((self.Ns, self.Nt))
        self.timepoints[self.masks["st"]] = timepoints_values

        d_values = np.concatenate(
            [hdxm.d_exp.to_numpy().flatten() for hdxm in self.hdxm_list]
        )
        self.d_exp = np.zeros((self.Ns, self.Np, self.Nt))
        self.d_exp[self.masks["spt"]] = d_values

        # Index array of of shape Ns x y where indices apply to dG return aligned residues for
        self.aligned_indices = None
        self.aligned_dataframes = None

    def __iter__(self):
        return self.hdxm_list.__iter__()

    def __getitem__(self, item: int) -> HDXMeasurement:
        return self.hdxm_list.__getitem__(item)

    def get(self, name: str) -> HDXMeasurement:
        """find a HDXMeasurement by name"""

        idx = self.names.index(name)
        return self[idx]

    @property
    def Ns(self) -> int:
        return len(self.hdxm_list)

    @property
    def Nr(self) -> int:
        return self.coverage.Nr

    @property
    def Np(self) -> int:
        return np.max([hdxm.Np for hdxm in self.hdxm_list])

    @property
    def Nt(self) -> int:
        return np.max([hdxm.Nt for hdxm in self.hdxm_list])

    @property
    def temperature(self) -> np.ndarray:
        return np.array([hdxm.temperature for hdxm in self.hdxm_list])

    @property
    def names(self) -> list[str]:
        return [hdxm.name for hdxm in self.hdxm_list]

    @property
    def rfu_residues(self) -> pd.DataFrame:
        """Relative fractional uptake per residue.

        Shape of the returned DataFrame is Nr (rows) x Ns*Nt (columns) and is multiindexed
        by columns (state, exposure, quantity)
        """
        rfu = pd.concat(
            [hdxm.rfu_residues for hdxm in self],
            keys=self.names,
            names=["state", "exposure"],
            axis=1,
        )
        columns = pd.MultiIndex.from_tuples(
            tuples=[(*tup, "rfu") for tup in rfu.columns],
            names=["state", "exposure", "quantity"],
        )

        rfu.columns = columns

        return rfu

    def guess_deltaG(self, rates_df: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """Obtain ΔG initial guesses from apparent H/D exchange rates.

        Args:
            rates_df: Pandas dataframe apparent exchange rates (units s^-1). Column names must
                correspond to HDX measurement names.
            **kwargs: Additional keyword arguments passed to :meth:`.HDXMeasurement.guess_deltaG`

        Returns:
            ΔG guess values (units kJ/mol)

        """

        guesses = [
            hdxm.guess_deltaG(rates_df[name], **kwargs)
            for hdxm, name in zip(self, self.names)
        ]
        deltaG = pd.concat(guesses, keys=self.names, axis=1)

        return deltaG

    # TODO alignment should be given as dict
    def add_alignment(self, alignment, first_r_numbers=None):
        """

        :param alignment: list
        :param first_r_numbers:
            default is [1, 1, ...] but specifiy here if alignments do not all start at residue 1
        :return:
        """
        dfs = [hdxm.coverage.protein.df for hdxm in self.hdxm_list]
        self.aligned_dataframes = align_dataframes(dfs, alignment, first_r_numbers)

        df = self.aligned_dataframes["r_number"]

        # Crop residue numbers to interval range
        df = df[
            ((self.coverage.interval[0] <= df) & (df < self.coverage.interval[1])).all(
                axis=1
            )
        ]
        df = (
            df - self.coverage.interval[0]
        )  # First residue in interval selected by index 0
        df.dropna(how="any", inplace=True)  # Remove non-aligned residues

        self.aligned_indices = df.to_numpy(dtype=int).T

    def get_tensors(
        self, dtype: Optional[torch.dtype] = None
    ) -> dict[str, torch.Tensor]:
        """Returns a dictionary of tensor variables for fitting HD kinetics.

        Tensor variables are (shape):
        Temperature (Ns x 1 x 1)
        X (Ns x Np x Nr)
        k_int (Ns x Nr)
        timepoints (Ns x 1 x Nt)
        d_exp (D) (Ns x Np x Nt)

        Returns:
            Dictionary with tensors

        """
        # todo create correct shapes as per table in docstring for all

        # TODO property?
        temperature = np.array([kf.temperature for kf in self.hdxm_list])

        X_values = np.concatenate(
            [hdxm.coverage.X.flatten() for hdxm in self.hdxm_list]
        )
        X = np.zeros((self.Ns, self.Np, self.Nr))
        X[self.masks["spr"]] = X_values

        k_int_values = np.concatenate(
            [hdxm.coverage["k_int"].to_numpy() for hdxm in self.hdxm_list]
        )
        k_int = np.zeros((self.Ns, self.Nr))
        k_int[self.masks["sr"]] = k_int_values

        dtype = dtype or cfg.TORCH_DTYPE
        device = cfg.TORCH_DEVICE

        tensors = {
            "temperature": torch.tensor(
                temperature, dtype=dtype, device=device
            ).reshape(self.Ns, 1, 1),
            "X": torch.tensor(X, dtype=dtype, device=device),
            "k_int": torch.tensor(k_int, dtype=dtype, device=device).reshape(
                self.Ns, self.Nr, 1
            ),
            "timepoints": torch.tensor(
                self.timepoints, dtype=dtype, device=device
            ).reshape(self.Ns, 1, self.Nt),
            "d_exp": torch.tensor(
                self.d_exp, dtype=dtype, device=device
            ),  # todo this is called uptake_corrected/D/uptake
        }

        return tensors

    @property
    def exchanges(self) -> np.ndarray:
        """Boolean mask ``True`` where there are residues which exchange

        Shape of the returned array is Ns x Np

        """
        values = np.concatenate(
            [hdxm.coverage["exchanges"].to_numpy() for hdxm in self.hdxm_list]
        )
        exchanges = np.zeros((self.Ns, self.Nr), dtype=bool)
        exchanges[self.masks["sr"]] = values

        return exchanges

    def to_file(
        self,
        file_path: os.PathLike,
        include_version: bool = True,
        include_metadata: bool = True,
        fmt: str = "csv",
        **kwargs: Any,
    ) -> None:
        """Write the data in this :class:`.HDXMeasurementSet` to file.

        Args:
            file_path: File path to create and write to.
            include_version: Set ``True`` to include PyHDX version and current time/date
            fmt: Formatting to use, options are 'csv' or 'pprint'
            include_metadata: If ``True``, the objects' metadata is included
            **kwargs: Optional additional keyword arguments passed to `df.to_csv`

        """

        dfs = []
        metadata = {}
        for hdxm in self.hdxm_list:
            metadata[hdxm.name] = (
                hdxm.metadata if include_metadata else include_metadata
            )
            dfs.append(hdxm.data)

        full_df = pd.concat(dfs, axis=1, keys=self.names)
        dataframe_to_file(
            file_path,
            full_df,
            include_version=include_version,
            include_metadata=metadata,
            fmt=fmt,
            **kwargs,
        )


# https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    (idx,) = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]  # Edit

    # Reshape the result into two columns
    idx.shape = (-1, 2)
    return idx


def hdx_intersection(
    hdx_list: list[HDXMeasurement], fields: Optional[list[str]] = None
):
    """
    Finds the intersection between peptides.

    Peptides are supplied as :class:`.HDXMeasurement` objects. After the intersection of
    peptides is found, new objects are returned where all peptides (coverage, exposure)
    between the measurements are identical.

    Optionally intersections by custom fields can be made.

    Args:
        hdx_list: Input list of :class:`.HDXMeasurement`
        fields: By which fields to take the intersections. Default is ['_start', '_end', 'exposure']

    Returns:
        hdx_out: Output list of :class:`.HDXMeasurement`
    """

    raise warnings.warn("'hdx_intersection' method is outdated", NotImplementedError)

    fields = fields or ["_start", "_end", "exposure"]

    full_arrays = [data_obj.full_data for data_obj in hdx_list]
    selected = array_intersection(full_arrays, fields=fields)

    hdx_out = [
        HDXMeasurement(data, **data_obj.metadata)
        for data, data_obj in zip(selected, hdx_list)
    ]
    return hdx_out


class PeptideUptakeModel(object):
    def __init__(self, sequence: list[str], temperature: float, pH: float) -> None:
        self.peptide = sequence[1:]  #
        self.temperature = temperature
        self.pH = pH
        padded_sequence = ["X"] + sequence + ["X"]
        k_int = k_int_from_sequence(padded_sequence, temperature, pH)
        self.k_int = k_int[2:-1]

    def eval_analytical(
        self, timepoints: np.ndarray, k_open: np.ndarray, k_close: np.ndarray
    ) -> np.ndarray:
        """Evaluate D-uptake for the given peptide at specified timepoints.


        Args:
            timepoints: Shape `(t,)` array with timepoints to sample.
            k_open: Shape `(k,)` array with opening rates (length equal to peptide length).
            k_close: Shape `(k,)` array with closing rates (length equal to peptide length).

        Returns:
            Shape (`t, k`) array with D-uptake values per amino acid per timepoint.
        """

        k_tot = k_open + k_close + self.k_int
        k_obs = 0.5 * (k_tot - np.sqrt((k_tot**2) - 4 * k_open * self.k_int))

        D_obs = 1 - np.exp(-k_obs[np.newaxis, :] * timepoints[:, np.newaxis])

        return D_obs

    def eval_single_numerical(
        self,
        aa_index: int,
        timepoints: np.ndarray,
        k_open: float,
        k_close: float,
        **solver_options: Any,
    ):

        k_tot = k_open + k_close
        y0 = np.array([k_close / k_tot, k_open / k_tot, 0])

        k_int = self.k_int[aa_index]
        if k_int == 0.0:
            return np.ones((len(timepoints), 3)) * y0

        method = solver_options.get("method", "LSODA")

        trs_rate_matrix = np.array(
            [[-k_open, k_close, 0], [k_open, -k_close - k_int, 0], [0, k_int, 0]]
        )

        jac = trs_rate_matrix if method != "LSODA" else None

        sol = solve_ivp(
            self.gradient_func,
            vectorized=True,
            t_span=(0, timepoints.max() * 1.001),
            y0=y0,
            jac=jac,
            t_eval=timepoints,
            args=[trs_rate_matrix],
            **solver_options,
        )

        return sol.y.T

    @staticmethod
    def gradient_func(t: Any, p: np.array, trs_matrix: np.array) -> np.ndarray:
        """
        calculates dp/dt given a transition state matrix and current populations p(at time t)

        :param p:
        :param t:
        :param trs_matrix: transition state matrix
        :return:
        """

        dpdt = trs_matrix @ p

        return dpdt

    def get_k_open(self, dG: npt.ArrayLike, k_close: npt.ArrayLike) -> npt.ArrayLike:
        return k_close / np.exp(dG / (R * self.temperature))

    def get_k_close(self, dG: npt.ArrayLike, k_open: npt.ArrayLike) -> npt.ArrayLike:
        return k_open * np.exp(dG / (R * self.temperature))

    def get_dG(self, k_open: npt.ArrayLike, k_close: npt.ArrayLike) -> npt.ArrayLike:
        return np.log(k_close / k_open) * (R * self.temperature)

    def __len__(self) -> int:
        return len(self.peptide)
