from __future__ import annotations

import os
import textwrap
import warnings
from numbers import Number
from typing import Optional, Any, Union, TYPE_CHECKING

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
from pyhdx.process import verify_sequence, parse_temperature, correct_d_uptake, apply_control
from pyhdx.support import reduce_inter, dataframe_intersection
from pyhdx.config import cfg

if TYPE_CHECKING:
    from hdxms_datasets import HDXDataSet


class Coverage:
    """
    Object describing layout and coverage of peptides and generating the corresponding matrices.
    Peptides should all belong to the same state and have the same exposure time.

    Args:
        data: DataFrame with input peptides
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
        self.protein = protein_df

        # matrix dimensions N_peptides N_residues, dtype for PyTorch compatibility
        _exchanges = self["exchanges"]  # Array only on covered part
        self.X = np.zeros((len(self.data), self.interval[1] - self.interval[0]), dtype=int)
        self.Z = np.zeros_like(self.X, dtype=float)
        for row, idx in enumerate(self.data.index):
            # start, end are already corrected for drop_first parameter
            start, end = self.data.loc[idx, "_start"], self.data.loc[idx, "_stop"]
            i0, i1 = self.r_number.get_loc(start), self.r_number.get_loc(end - 1)
            # i0, i1 = np.searchsorted(self.r_number, (entry['start'], entry['end']))
            self.X[row][i0 : i1 + 1] = 1
            self.Z[row][i0 : i1 + 1] = _exchanges.iloc[i0 : i1 + 1]
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

    def apply_interval(self, array_or_series: Union[np.ndarray, pd.Series]) -> pd.Series:
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
        """Pandas index of residue numbers corresponding to the part of the protein covered by peptides."""
        return self.r_number

    @property
    def block_length(self) -> np.ndarray:
        """Lengths of unique blocks of residues in the peptides map, along the `r_number` axis"""

        # indices are start and stop values of blocks
        indices = np.sort(np.concatenate([self.data["_start"], self.data["_stop"]]))
        # indices of insertion into r_number vector gives us blocks with taking prolines into account.
        diffs = np.diff(np.searchsorted(self.r_number, indices))

        block_length = diffs[diffs != 0]
        return block_length

    @property
    def X_norm(self) -> np.ndarray:
        """`X` coefficient matrix normalized column-wise."""
        return self.X / np.sum(self.X, axis=0)[np.newaxis, :]

    @property
    def Z_norm(self) -> np.ndarray:
        """`Z` Coefficient matrix normalized column-wise."""
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


class HDXMeasurement:
    """Main HDX data object.

    This object has peptide data of a single state and with multiple timepoints.
    Timepoint data is split into [`HDXTimepoint`][models.HDXTimepoint] objects for
    each timepoint. Supplied data is made 'uniform' such that all timepoints have the same peptides.

    Args:
        data: Dataframe with all peptides belonging to a single state.
        **metadata: Dictionary of optional metadata. By default, holds the `temperature` and `pH` parameters.

    Attributes:
        coverage: Coverage object describing peptide layout.
        data: DataFrame with all peptides, taking only peptides present in all timepoints.
        peptides: List of `HDXTimepoint` objects, one per exposure.
        state: Protein state label for this HDX measurement.
        timepoints: Deuterium exposure times.
    """

    def __init__(self, data: pd.DataFrame, **metadata: Any):
        self.metadata = metadata
        assert len(data["state"].unique()) == 1
        self.state: str = str(data["state"].iloc[0])
        self.timepoints: np.ndarray = np.sort(np.unique(data["exposure"]))

        # todo sort happens twice now
        data = data.sort_values(["start", "stop", "sequence", "exposure"])

        # Obtain the intersection of peptides per timepoint
        df_list = [(data[data["exposure"] == exposure]) for exposure in self.timepoints]

        intersected_data = dataframe_intersection(df_list, by=["start", "stop"])

        cov_kwargs = {kwarg: metadata.get(kwarg) for kwarg in ["c_term", "n_term", "sequence"]}
        self.peptides: list[HDXTimepoint] = [
            HDXTimepoint(df, **cov_kwargs) for df in intersected_data
        ]

        # Create coverage object from the first time point (as all are now equal)
        self.coverage: Coverage = Coverage(intersected_data[0], **cov_kwargs)

        if self.temperature and self.pH:
            # list(self.protein["sequence"])
            k_int_array = k_int_from_sequence(
                self.coverage.protein["sequence"], self.temperature, self.pH
            )

            # k_int = self.coverage.protein.get_k_int(self.temperature, self.pH)
            self.coverage.protein["k_int"] = k_int_array

        self.data: pd.DataFrame = pd.concat(
            intersected_data, axis=0, ignore_index=True
        ).sort_values(["start", "stop", "sequence", "exposure"])
        self.data["peptide_id"] = self.data.index % self.Np
        self.data.index.name = (
            "peptide_index"  # index is original index which continues along exposures
        )
        self.data_wide = (
            self.data.pivot(index="peptide_id", columns=["exposure"])
            .reorder_levels([1, 0], axis=1)
            .sort_index(axis=1, level=0, sort_remaining=False)
        )

    @classmethod
    def from_dataset(cls, dataset: HDXDataSet, state: str | int, **metadata) -> HDXMeasurement:
        """Create an HDXMeasurement object from a HDXDataSet object.

        Args:
            dataset: HDXDataSet object
            state: State label or index for measurement in the dataset

        Returns:
            HDXMeasurement object.

        """

        state = dataset.states[state] if isinstance(state, int) else state
        peptide_spec = dataset.hdx_spec["states"][state]["peptides"]

        peptides = dataset.load_peptides(state, "experiment")
        if "FD_control" not in peptide_spec:
            raise ValueError("Dataset does not contain a FD_control state")
        fd_peptides = dataset.load_peptides(state, "FD_control")
        nd_peptides = (
            dataset.load_peptides(state, "ND_control") if "ND_control" in peptide_spec else None
        )

        # take globally defined metadata and update with state specific metadata
        spec_metadata = dataset.hdx_spec.get("metadata", {})
        spec_metadata.update(dataset.hdx_spec["states"][state]["metadata"])

        metadata = {**spec_metadata, **metadata}

        peptides = apply_control(peptides, fd_peptides, nd_peptides)
        peptides = correct_d_uptake(
            peptides,
            drop_first=cfg.analysis.drop_first,
            d_percentage=metadata.get("d_percentage", 100.0),
        )

        return HDXMeasurement(peptides, name=state, **metadata)

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
        """HDX Measurement name."""
        return self.metadata.get("name", self.state)

    @property
    def temperature(self) -> Optional[float]:
        """Temperature of the H/D exchange reaction (K)."""
        temperature = self.metadata.get("temperature")
        if isinstance(temperature, (Number, type(None))):
            return temperature
        elif isinstance(temperature, dict):
            return parse_temperature(**temperature)

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

    def __iter__(self):
        return self.peptides.__iter__()

    def __getitem__(self, item):
        return self.peptides.__getitem__(item)

    @property
    def rfu_residues(self) -> pd.DataFrame:
        """Relative fractional uptake per residue.

        Shape of the returned DataFrame is `(Nr, Nt)`.
        """
        df = pd.concat([v.rfu_residues for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"

        return df

    @property
    def rfu_residues_sd(self) -> pd.DataFrame:
        """Standard deviations of relative fractional uptake per residue.

        Shape of the returned DataFrame is `(Nr, Nt)`.
        """

        df = pd.concat([v.rfu_residues_sd for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"

        return df

    @property
    def rfu_peptides(self) -> pd.DataFrame:
        """Relative fractional uptake per peptide.

        Shape of the returned DataFrame is `(Np, Nt)`.
        """
        df = pd.concat([v.rfu_peptides for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"
        return df

    @property
    def d_exp(self) -> pd.DataFrame:
        """D-uptake values (corrected for back-exchange).

        Shape of the returned DataFrame is `(Np, Nt)`.
        """
        df = pd.concat([v.d_exp for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"
        return df

    # todo check shapes of k_int and timepoints, compared to their shapes in hdxmeasurementset
    def get_tensors(
        self, exchanges: bool = False, dtype: Optional[torch.dtype] = cfg.TORCH_DTYPE
    ) -> dict[str, torch.Tensor]:
        """Returns a dictionary of tensor variables for fitting HD kinetics.

        Args:
            exchanges: If `True` only returns tensor data describing residues which exchange
                (ie have peptides and are not prolines).
            dtype: Optional Torch data type. Use torch.float32 for faster fitting of large data
                sets, possibly at the expense of accuracy.

        Returns:
            Dictionary with tensors.

        Note: Tensor output and shapes:
            * temperature `(1, 1)`
            * X `(Np, Nr)`
            * k_int `(Nr, 1)`
            * timepoints `(1, Nt)`
            * d_exp `(Np, Nt)`
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
            "temperature": torch.tensor([self.temperature], dtype=dtype, device=device).unsqueeze(
                -1
            ),
            "X": torch.tensor(self.coverage.X[:, bools], dtype=dtype, device=device),
            "k_int": torch.tensor(
                self.coverage["k_int"].to_numpy()[bools], dtype=dtype, device=device
            ).unsqueeze(-1),
            "timepoints": torch.tensor(self.timepoints, dtype=dtype, device=device).unsqueeze(0),
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
            rates: Apparent exchange rates (units s^-1^). Series index is protein residue number.
            correct_c_term: If `True`, sets the guess value of the c-terminal residue to the
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

        p_guess.clip(0.0, None, inplace=True)  # Some initial guesses might have negative PF values
        with np.errstate(divide="ignore"):
            deltaG = np.log(p_guess) * constants.R * self.temperature

        deltaG.replace([np.inf, -np.inf], np.nan, inplace=True)

        c_term = self.coverage.protein.index.max()
        if correct_c_term and c_term in deltaG.index:
            deltaG.loc[c_term] = deltaG.loc[c_term - 1]

        return deltaG

    def to_file(
        self,
        file_path: os.PathLike,
        include_version: bool = True,
        include_metadata: bool = True,
        fmt: str = "csv",
        **kwargs: Any,
    ) -> None:
        """Write the data in this [HDXMeasurement][models.HDXMeasurement] to file.

        Args:
            file_path: File path to create and write to.
            include_version: Set to `True` to include PyHDX version and current time/date
            fmt: Formatting to use, options are 'csv' or 'pprint'
            include_metadata: If `True`, the objects' metadata is included
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
    """Class with subset of peptides corresponding to only one state and exposure.

    Args:
        data: Dataframe with input data.
        **kwargs: Additional keyword arguments passed to [Coverage][models.Coverage].

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

        return self.propagate_errors("rfu_sd")

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
        # TODO this should be functional; remove coverage object as a whole
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


class CoverageSet:
    """Coverage object for multiple [HDXMeasurement][models.HDXMeasurement] objects.

    This objects finds the minimal interval of residue numbers which fit all :class:`.HDXMeasurement`s


    Args:
        hdxm_list: List of input HDXMeasurment objects.

    """

    # todo perhaps this object should have X
    def __init__(self, hdxm_list: list[HDXMeasurement]):
        self.hdxm_list = hdxm_list

        # todo create Coverage object for the 3d case
        intervals = np.array([hdxm_list.coverage.interval for hdxm_list in self.hdxm_list])
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

    def get_masks(self) -> dict[str, np.ndarray]:
        """Get boolean masks along the different data dimensions which are `True` at elements
            which have measured data.

        The masks can be used to assign values to the quantity tensors (k_int, X, D_exp, ect) used
        for calculating D-uptake from ΔG.

        Returns:
            Dictionary of boolean masks spanning the various data dimensions.

        Note: Returned masks and shapes:
            * sr_mask: `(Ns, Nr)`
            * st_mask: `(Ns, Nt)`
            * spt_mask: `(Ns, Np, Nr)`
            * spt_mask: `(Ns, Np, Nt)`

        """
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


class HDXMeasurementSet:
    """
    Set of multiple [HDXMeasurement][models.HDXMeasurement] objects.

    Args:
        hdxm_list: Input list of HDX measurements.

    Attributes:
        coverage: Coverage object for the set of measurements.
        d_exp: Array with measured D-uptake values, padded with zeros. Shape  is `(Ns, Np, Nt)`.
        timepoints: Array with timepoints, padded with zeros in case of samples with
            unequal number of timepoints. Shape is `(Ns, Nt)`.

    """

    def __init__(self, hdxm_list: list[HDXMeasurement]) -> None:
        self.hdxm_list = hdxm_list

        self.coverage: CoverageSet = CoverageSet(hdxm_list)
        self.masks = self.coverage.get_masks()

        timepoints_values = np.concatenate([hdxm.timepoints for hdxm in self.hdxm_list])
        self.timepoints: np.ndarray = np.zeros((self.Ns, self.Nt))
        self.timepoints[self.masks["st"]] = timepoints_values

        d_values = np.concatenate([hdxm.d_exp.to_numpy().flatten() for hdxm in self.hdxm_list])
        self.d_exp: np.ndarray = np.zeros((self.Ns, self.Np, self.Nt))
        self.d_exp[self.masks["spt"]] = d_values

        # Index array of shape Ns x y where indices apply to dG return aligned residues for
        self.aligned_indices = None
        self.aligned_dataframes = None

    def __iter__(self):
        return self.hdxm_list.__iter__()

    def __getitem__(self, item: int) -> HDXMeasurement:
        return self.hdxm_list.__getitem__(item)

    @classmethod
    def from_dataset(self, dataset: HDXDataSet, **metadata) -> HDXMeasurementSet:
        hdxm_list = [
            HDXMeasurement.from_dataset(dataset, state, **metadata) for state in dataset.states
        ]

        return HDXMeasurementSet(hdxm_list)

    def get(self, name: str) -> HDXMeasurement:
        """
        Get HDXMeasurement object by name.

        Args:
            name: Name of the HDXMeasurement object.

        Returns:
            The HDXMeasurement object
        """

        idx = self.names.index(name)
        return self[idx]

    @property
    def Ns(self) -> int:
        """Number of samples"""
        return len(self.hdxm_list)

    @property
    def Nr(self) -> int:
        """Number of residues"""
        return self.coverage.Nr

    @property
    def Np(self) -> int:
        """Number of peptides"""
        return np.max([hdxm.Np for hdxm in self.hdxm_list])

    @property
    def Nt(self) -> int:
        """Number of timepoints"""
        return np.max([hdxm.Nt for hdxm in self.hdxm_list])

    @property
    def temperature(self) -> np.ndarray:
        """Array of temperature values for each measurement"""
        return np.array([hdxm.temperature for hdxm in self.hdxm_list])

    @property
    def names(self) -> list[str]:
        """List of names of the measurement"""
        return [hdxm.name for hdxm in self.hdxm_list]

    @property
    def rfu_residues(self) -> pd.DataFrame:
        """Relative fractional uptake per residue.

        Returned DataFrame has shape `(Nr, Ns*Nt)`, which is multiindex by columns (state, exposure, quantity).

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

    def guess_deltaG(self, rates_df: pd.DataFrame, correct_c_term: bool = True) -> pd.DataFrame:
        """Obtain ΔG initial guesses from apparent H/D exchange rates.

        Args:
            rates_df: Pandas dataframe apparent exchange rates (units s^-1). Column names must
                correspond to HDX measurement names.
            correct_c_term: If `True`, sets the guess value of the c-terminal residue to the
                value of the residue preceding it.

        Returns:
            ΔG guess values (units J/mol).
        """

        guesses = [
            hdxm.guess_deltaG(rates_df[name], correct_c_term=correct_c_term)
            for hdxm, name in zip(self, self.names)
        ]
        deltaG = pd.concat(guesses, keys=self.names, axis=1)

        return deltaG

    # TODO alignment should be given as dict
    def add_alignment(self, alignment, first_r_numbers=None) -> None:
        """
        Args:
            alignment: FASTA alignments.
            first_r_numbers: default is [1, 1, ...] but specifiy here if alignments do not all start at residue 1

        """

        dfs = [hdxm.coverage.protein for hdxm in self.hdxm_list]
        self.aligned_dataframes = align_dataframes(dfs, alignment, first_r_numbers)

        df = self.aligned_dataframes["r_number"]

        # Crop residue numbers to interval range
        df = df[((self.coverage.interval[0] <= df) & (df < self.coverage.interval[1])).all(axis=1)]
        df = df - self.coverage.interval[0]  # First residue in interval selected by index 0
        df.dropna(how="any", inplace=True)  # Remove non-aligned residues

        self.aligned_indices = df.to_numpy(dtype=int).T

    def get_tensors(self, dtype: Optional[torch.dtype] = None) -> dict[str, torch.Tensor]:
        """Returns a dictionary of tensor variables for fitting HD kinetics.

        Args:
            dtype: Optional Torch data type. Use torch.float32 for faster fitting of large data
                sets, possibly at the expense of accuracy.

        Returns:
            Dictionary with tensors.

        Note: Tensor output and shapes:
            * temperature `(Ns, 1, 1)`
            * X `(Ns, Np, Nr)`
            * k_int `(Ns, Nr, 1)`
            * timepoints `(Ns, 1, Nt)`
            * d_exp `(Ns, Np, Nt)`

        """
        # todo create correct shapes as per table in docstring for all

        # TODO property?
        temperature = np.array([kf.temperature for kf in self.hdxm_list])

        X_values = np.concatenate([hdxm.coverage.X.flatten() for hdxm in self.hdxm_list])
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
            "temperature": torch.tensor(temperature, dtype=dtype, device=device).reshape(
                self.Ns, 1, 1
            ),
            "X": torch.tensor(X, dtype=dtype, device=device),
            "k_int": torch.tensor(k_int, dtype=dtype, device=device).reshape(self.Ns, self.Nr, 1),
            "timepoints": torch.tensor(self.timepoints, dtype=dtype, device=device).reshape(
                self.Ns, 1, self.Nt
            ),
            "d_exp": torch.tensor(
                self.d_exp, dtype=dtype, device=device
            ),  # todo this is called uptake_corrected/D/uptake
        }

        return tensors

    @property
    def exchanges(self) -> np.ndarray:
        """Boolean mask for residues which exchange (shape `(Ns, Np)`)"""
        values = np.concatenate([hdxm.coverage["exchanges"].to_numpy() for hdxm in self.hdxm_list])
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
            metadata[hdxm.name] = hdxm.metadata if include_metadata else include_metadata
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


class PeptideUptakeModel:
    """Model D-uptake in a single peptide.

    Args:
        sequence: FASTA sequence as list of strings.
        temperature: Temperature of the H/D exchange reaction in Kelvin.
        pH: pH of the H/D exchange reaction.

    Attributes:
        k_int: Array of intrinsic exchanges rates


    """

    def __init__(self, sequence: list[str], temperature: float, pH: float) -> None:
        self.peptide = sequence[1:]  #
        self.temperature = temperature
        self.pH = pH
        padded_sequence = ["X"] + sequence + ["X"]
        k_int = k_int_from_sequence(padded_sequence, temperature, pH)
        self.k_int: np.array = k_int[2:-1]

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
