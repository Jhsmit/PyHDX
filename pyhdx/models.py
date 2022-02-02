import textwrap
import warnings
from functools import reduce, partial

import numpy as np
import pandas as pd
import torch
from hdxrate import k_int_from_sequence
from numpy.lib.recfunctions import append_fields
from scipy import constants

import pyhdx
from pyhdx.alignment import align_dataframes
from pyhdx.fileIO import dataframe_to_file
from pyhdx.support import reduce_inter, fields_view
from pyhdx.config import cfg


def protein_wrapper(func, *args, **kwargs):
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


class PeptideMasterTable(object):
    """
    Main peptide input object.

    The input `~pandas.DataFrame`pandas DataFrame `data` must have the following
    entires for each peptide:

    start: Residue number of the first amino acid in the peptide
    end: Residue number of the last amino acid in the peptide (inclusive)
    sequence: Amino acid sequence of the peptide (one letter code)
    exposure: Typically the time the sample was exposed to a deuterated solution. This can correspond to other times if
        the kinetics of the experiment are set up differently
    state: String describing to which state (experimental conditions) the peptide belongs
    uptake: Number of deuteriums the peptide has taken up

    The following fields are added to the `data` array upon initialization:

    _start: Unmodified copy of initial start field
    _end: Unmodified copy of initial end field
    _sequence: Unmodified copy of initial sequence
    ex_residues: Number of residues that undergo deuterium exchange. This number is calculated using the `drop_first` and
        `ignore_prolines` parameters

    N-terminal residues which are removed because they are either within `drop_first` or they are N-terminal prolines are
    marked with 'x' in the `sequence` field. Prolines which are removed because they are in the middle of a peptide are
    marked with a lower case 'p' in the sequence field.

    The field `scores` is used in calculating exchange rates and can be set by either the `set_backexchange` or
    `set_control` methods.


    Parameters
    ----------
    data : :class:`~pandas.DataFrame`
        Pandas DataFrame with peptide entries.
    drop_first : :obj:`int`
        Number of N-terminal amino acids to ignore. Default is 1.
    d_percentage : :obj:`float`
        Percentage of deuterium in the labelling solution.
    ignore_prolines : :obj:`bool`
        Boolean to toggle ignoring of proline residues. When True these residues are treated as if they're not present
        in the protein.
    sort : :obj:`bool`
        Set to ``True`` to sort the input. Sort order is 'start', 'end', 'sequence', 'exposure', 'state'.
    remove_nan : :obj:`bool`
        Set to ``True`` to remove NaN entries in uptake

    """

    def __init__(
        self,
        data,
        drop_first=1,
        ignore_prolines=True,
        d_percentage=100.0,
        sort=True,
        remove_nan=True,
    ):
        assert np.all(
            data["start"] < data["end"]
        ), "All `start` entries must be smaller than their `end` entries"
        assert (
            0 <= d_percentage <= 100.0
        ), "Deuteration percentage must be between 0 and 100"
        d_percentage /= 100.0

        self.data = data.copy().reset_index(drop=True)
        self.data.index.name = "peptide_index"

        if remove_nan:
            self.data = self.data.dropna(subset=["uptake"])
        if sort:
            self.data = self.data.sort_values(
                ["start", "end", "sequence", "state", "exposure"]
            )

        for col in ["start", "end", "sequence"]:
            target = "_" + col
            if target in self.data:
                continue
            else:
                self.data[target] = self.data[col]

        # Convert sequence to upper case if not so already
        self.data["sequence"] = self.data["sequence"].str.upper()
        # Mark ignored prolines with lower case letters
        if ignore_prolines:
            self.data["sequence"] = [s.replace("P", "p") for s in self.data["sequence"]]

        # Find the total number of n terminal / c_terminal residues to remove
        # Todo: edge cases such as pure prolines or overlap between c terminal prolines and drop_first section (issue 32)
        n_term = np.array(
            [
                len(seq) - len(seq[drop_first:].lstrip("p"))
                for seq in self.data["sequence"]
            ]
        )
        c_term = np.array(
            [len(seq) - len(seq.rstrip("p")) for seq in self.data["sequence"]]
        )

        # Mark removed n terminal residues with lower case x
        self.data["sequence"] = [
            "x" * nt + s[nt:] for nt, s in zip(n_term, self.data["sequence"])
        ]
        self.data["start"] += n_term
        self.data["end"] -= c_term

        ex_residues = (
            np.array(
                [len(s) - s.count("x") - s.count("p") for s in self.data["sequence"]]
            )
            * d_percentage
        )
        if "ex_residues" not in self.data:
            self.data["ex_residues"] = ex_residues

    def __len__(self):
        return self.data.shape[0]

    def get_state(self, state):
        """
        Returns entries in the table with state 'state'
        Rows with NaN entries for 'uptake_corrected' are removed

        Parameters
        ----------
        state : :obj:`str`


        Returns
        -------

        """

        if not isinstance(state, str):
            raise TypeError(f"State must be type `str`, got {type(state)}")
        data = self.data.query(f'state == "{state}"').copy()
        if "uptake_corrected" in data.columns:
            data.dropna(subset=["uptake_corrected"], inplace=True)
        if len(data) == 0:
            raise ValueError(
                f"No data found for state {state!r}, options are: {', '.join(self.data['state'].unique())}"
            )

        return data

    def set_backexchange(self, back_exchange):
        """
        Sets the normalized percentage of uptake through a fixed backexchange value for all peptides.

        Parameters
        ----------
        back_exchange :  :obj:`float`
            Percentage of back exchange

        """

        back_exchange /= 100
        rfu = self.data["uptake"] / ((1 - back_exchange) * self.data["ex_residues"])

        uptake_corrected = self.data["uptake"] / (1 - back_exchange)

        self.data = append_fields(
            self.data,
            ["rfu", "uptake_corrected"],
            data=[rfu, uptake_corrected],
            usemask=False,
        )

    def set_control(self, control_1, control_0=None):
        """
        Apply a control dataset to this object. The column 'RFU' is added to the object by normalizing its uptake
        value with respect to the control uptake value to one.
        Optionally, ``control_zero`` can be specified which is a dataset whose uptake value will be used to zero
        the uptake.

        Nonmatching peptides are set to NaN

        #todo insert math

        Parameters
        ----------
        control_1 : :obj:`tuple`
            tuple with (`state`, `exposure`) for peptides to use for normalization (FD control)
        control_0 : :obj:`tuple`, optional
            tuple with (`state`, `exposure`) for peptides to use for zeroing uptake values (ND control)

        """

        try:
            fd_df = self.get_data(*control_1)[["_start", "_end", "uptake"]].set_index(
                ["_start", "_end"], verify_integrity=True
            )
        except ValueError as e:
            raise ValueError("FD control has duplicate entries") from e

        if fd_df.size == 0:
            raise ValueError(
                f"No matching peptides with state {control_1[0]} and exposure {control_1[1]}"
            )

        try:
            if control_0 is None:
                nd_df = (
                    self.get_data(*control_1)
                    .copy()[["_start", "_end", "uptake"]]
                    .set_index(["_start", "_end"], verify_integrity=True)
                )
                nd_df["uptake"] = 0

            else:
                nd_df = self.get_data(*control_0)[
                    ["_start", "_end", "uptake"]
                ].set_index(["_start", "_end"], verify_integrity=True)
                if nd_df.size == 0:
                    raise ValueError(
                        f"No matching peptides with state {control_0[0]} and exposure {control_0[1]}"
                    )
        except ValueError as e:
            raise ValueError("ND control has duplicate entries") from e

        self.data.set_index(["_start", "_end"], append=True, inplace=True)
        self.data.reset_index(level=0, inplace=True)

        self.data["rfu"] = (self.data["uptake"] - nd_df["uptake"]) / (
            fd_df["uptake"] - nd_df["uptake"]
        )
        self.data["uptake_corrected"] = self.data["rfu"] * self.data["ex_residues"]

        self.data = self.data.set_index("peptide_index", append=True).reset_index(
            level=[0, 1]
        )

    def select(self, **kwargs):
        """
        Select data based on column values.

        Parameters
        ----------
        kwargs: :obj:`dict`
            Column name, value pairs to select

        Returns
        -------
        output_data : :class:`~pandas.DataFrame`
            DataFrame with selected peptides

        """
        masks = [self.data[k] == v for k, v in kwargs.items()]
        m = np.logical_and.reduce(masks)

        return self.data[m]

    def get_data(self, state, exposure):
        """
        Get all peptides matching `state` and `exposure`.

        Parameters
        ----------
        state : :obj:`str`
            Measurement state
        exposure : :obj:`float`
            Measurement exposure time

        Returns
        -------
        output_data : :class:`~pandas.DataFrame`
            DataFrame with selected peptides
        """

        return self.select(state=state, exposure=exposure)

    @property
    def states(self):
        """:class:`~numpy.ndarray` Array with unique states"""
        return np.unique(self.data["state"])

    @property
    def exposures(self):
        """:class:`~numpy.ndarray` Array with unique exposures"""
        return np.unique(self.data["exposure"])


class Coverage(object):
    """
    Object describing layout and coverage of peptides and generating the corresponding matrices. Peptides should all
    belong to the same state and have the same exposure time.

    Parameters
    ----------
    data : :class:`~pandas.DataFrame`
        DataFrame with input peptides
    c_term : :obj:`int`
        Residue index number of the C-terminal residue (where first residue in index number 1)
    n_term : :obj:`int`
        Residue index of the N-terminal residue. Default value is 1, can be negative to accomodate for N-terminal
        purification tags
    sequence : :obj:`str`
        Amino acid sequence of the protein in one-letter FASTA encoding. Optional, if not specified the amino acid sequence
        from the peptide data is used to (partially) reconstruct the sequence. Supplied amino acid sequence must be
        compatible with sequence information in the peptides.

    Attributes
    ----------

    X : :class:`~numpy.ndarray`
        N x M matrix where N is the number of peptides and M equal to `prot_len`.
        Values are 1/(ex_residues) where there is coverage.
    Z : :class:`~numpy.ndarray`
        N x M matrix where N is the number of peptides and M equal to `prot_len`.
        Values are 1/(ex_residues) where there is coverage,
        #todo account for prolines: so that rows sum to 1 is currently not true

    """

    def __init__(self, data, n_term=1, c_term=None, sequence=""):
        assert len(np.unique(data["exposure"])) == 1, "Exposure entries are not unique"
        assert len(np.unique(data["state"])) == 1, "State entries are not unique"
        self.data = data.sort_values(["start", "end"], axis=0)
        self.data.index.name = "peptide_id"  # todo check these are the same as parent object peptide_id (todo make wide instead of instersection)

        start = self.data["_start"].min()
        end = self.data["_end"].max()

        if n_term:
            start = min(start, n_term)
        if sequence and c_term is None:
            c_term = len(sequence) + n_term - 1
        if c_term:
            if c_term + 1 < end:
                raise ValueError(
                    "HDX data extends beyond c_term number, check 'sequence' or 'c_term'"
                )
            end = c_term + 1  # c_term is inclusive, therefore plus one

        r_number = pd.RangeIndex(
            start, end, name="r_number"
        )  # r_number spanning the full protein range, not just the covered range
        # Full sequence
        _seq = pd.Series(index=r_number, dtype="U").fillna("X")  # Full sequence
        # Sequence with lower case letters for no coverage due to n_terminal residues or prolines
        seq = pd.Series(index=r_number, dtype="U").fillna("X")
        for idx in self.data.index[::-1]:
            start, end = self.data.loc[idx, "_start"], self.data.loc[idx, "_end"]

            _seq.loc[start : end - 1] = list(self.data.loc[idx, "_sequence"])
            seq.loc[start : end - 1] = list(
                self.data.loc[idx, "sequence"]
            )  # = list(d['sequence'])

        if sequence:
            for r, s1, s2 in zip(r_number, sequence, _seq):
                if s2 != "X" and s1 != s2:
                    raise ValueError(
                        f"Mismatch in supplied sequence and peptides sequence at residue {r}, expected '{s2}', got '{s1}'"
                    )
            if len(sequence) != len(_seq):
                raise ValueError(
                    "Invalid length of supplied sequence. Please check 'n_term' and 'c_term' parameters"
                )
            _seq = list(sequence)

        # todo check if this is always correctly determined (n terminal residues usw)
        exchanges = [
            s.isupper() and (s != "X") for s in seq
        ]  # Boolean array True if residue exchanges, full length
        coverage = seq != "X"  # Boolean array for coverage
        protein_df = pd.DataFrame(
            {"sequence": _seq, "coverage": coverage, "exchanges": exchanges},
            index=r_number,
        )

        # Inclusive, exclusive interval of peptides coverage across the whole protein
        self.interval = (np.min(self.data["start"]), np.max(self.data["end"]))
        self.protein = Protein(protein_df, index="r_number")

        # matrix dimensions N_peptides N_residues, dtype for TF compatibility
        _exchanges = self["exchanges"]  # Array only on covered part
        self.X = np.zeros(
            (len(self.data), self.interval[1] - self.interval[0]), dtype=int
        )
        self.Z = np.zeros_like(self.X, dtype=float)
        for row, idx in enumerate(self.data.index):
            start, end = self.data.loc[idx, "start"], self.data.loc[idx, "end"]
            i0, i1 = self.r_number.get_loc(start), self.r_number.get_loc(end - 1)
            # i0, i1 = np.searchsorted(self.r_number, (entry['start'], entry['end']))
            self.X[row][i0 : i1 + 1] = 1
            self.Z[row][i0 : i1 + 1] = _exchanges[i0 : i1 + 1]
        self.Z = self.Z / self.data["ex_residues"].to_numpy()[:, np.newaxis]

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        pd_series = self.protein[item]
        return self.apply_interval(pd_series)

    def apply_interval(self, array_or_series):
        """
        Given a Numpy array or Pandas series with a length equal to the full protein, returns the section of the array equal to the covered
        region. Returned series length is equal to number of columns in the X matrix

        Parameters
        ----------
        np.narray or pd.series

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
    def percent_coverage(self):
        """:obj:`float`: Percentage of residues covered by peptides"""
        return 100 * np.mean(self.protein["coverage"])

    @property
    def redundancy(self):
        """:obj:`float`: Average redundancy of peptides in regions with at least 1 peptide"""
        x_coverage = self.X[:, self["coverage"]]
        return np.mean(np.sum(x_coverage, axis=0))

    @property
    def Np(self):
        """:obj:`int`: Number of peptides."""
        return self.X.shape[0]

    @property
    def Nr(self):
        """:obj:`int`: Total number of residues spanned by the peptides."""

        return self.X.shape[1]

    @property
    def r_number(self):
        """:class:`~pandas.RangeIndex`: Pandas index numbers corresponding to the part of the protein covered by peptides"""

        return pd.RangeIndex(self.interval[0], self.interval[1], name="r_number")

    @property
    def index(self):
        """:class:`~pandas.RangeIndex`: Pandas index numbers corresponding to the part of the protein covered by peptides"""
        return self.r_number

    @property
    def block_length(self):
        """:class:`~numpy.ndarary`: Lengths of unique blocks of residues in the peptides map,
        along the `r_number` axis"""

        # indices are start and stop values of blocks
        indices = np.sort(np.concatenate([self.data["start"], self.data["end"]]))
        # indices of insertion into r_number vector gives us blocks with taking prolines into account.
        diffs = np.diff(np.searchsorted(self.r_number, indices))

        block_length = diffs[diffs != 0]
        return block_length

    @property
    def X_norm(self):
        """:class:`~numpy.ndarray`: `X` coefficient matrix normalized column wise."""
        return self.X / np.sum(self.X, axis=0)[np.newaxis, :]

    @property
    def Z_norm(self):
        """:class:`~numpy.ndarray`: `Z` coefficient matrix normalized column wise."""
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            z_norm = self.Z / np.sum(self.Z, axis=0)[np.newaxis, :]
        return z_norm

    def get_sections(self, gap_size=-1):
        """Get the intervals of independent sections of coverage.

        Intervals are inclusive, exclusive.
        Gaps are defined with `gap_size`, adjacent peptides with distances bigger than this value are considered not to
        overlap. Set to -1 to treat touching peptides as belonging to the same section.

        Parameters
        ----------
        gap_size : :obj:`int`
            The size which defines a gap

        """
        intervals = [(s, e) for s, e in zip(self.data["start"], self.data["end"])]
        sections = reduce_inter(intervals, gap_size=gap_size)

        return sections

    def __eq__(self, other):
        """Coverage objects are considered equal if both objects fully match between their start, end and sequence fields"""
        assert isinstance(other, Coverage), "Other must be an instance of Coverage"
        return (
            len(self.data) == len(other.data)
            and np.all(self.data["start"] == other.data["start"])
            and np.all(self.data["end"] == other.data["end"])
            and np.all(self.data["sequence"] == other.data["sequence"])
        )


class HDXMeasurement(object):
    """
    Main HDX data object. This object has peptide data of a single state but with multiple timepoints.

    Timepoint data is split into :class:`~pyhdx.models.PeptideMeasurements` objects for each timepoint
    Supplied data is made 'uniform' such that all timepoints have the same peptides

    Parameters
    ----------
    data : :class:`~pandas.DataFrame`
        Pandas dataframe with all peptides belonging to a single state.
    **metadata
        Dictionary of optional metadata. By default, holds the `temperature` and `pH` parameters.


    Attributes
    ----------
    data: :class:`~pandas.DataFrame`
        Pandas dataframe with all peptides
    state : :obj:`str`
        State of the HDX measurement
    timepoints : :class:`~numpy.ndarray`
        Array with exposure times (sorted)
    peptides : :obj:`list`
        List of :class:`~pyhdx.models.PeptideMeasurements`, one list element per timepoint.
    coverage : :class:`~pyhdx.models.Coverage`
        Coverage object describing peptide layout.

    """

    def __init__(self, data, **metadata):
        self.metadata = metadata
        assert len(data["state"].unique()) == 1
        self.state = str(data["state"].iloc[0])
        self.timepoints = np.sort(np.unique(data["exposure"]))

        # Obtain the intersection of peptides per timepoint
        data_list = [
            (data[data["exposure"] == exposure]).set_index(["_start", "_end"])
            for exposure in self.timepoints
        ]
        index_intersection = reduce(pd.Index.intersection, [d.index for d in data_list])
        intersected_data = [
            df.loc[index_intersection].reset_index() for df in data_list
        ]

        cov_kwargs = {
            kwarg: metadata.get(kwarg, default)
            for kwarg, default in zip(["c_term", "n_term", "sequence"], [None, 1, ""])
        }

        self.peptides = [HDXTimepoint(df, **cov_kwargs) for df in intersected_data]

        # Create coverage object from the first time point (as all are now equal)
        self.coverage = Coverage(intersected_data[0], **cov_kwargs)

        if self.temperature and self.pH:
            k_int = self.coverage.protein.get_k_int(self.temperature, self.pH)
            self.coverage.protein["k_int"] = k_int

        self.data = pd.concat(intersected_data, axis=0, ignore_index=True).sort_values(
            ["start", "end", "sequence", "exposure"]
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

    def __str__(self):
        """
        String representation of HDX measurement object.

        Returns
        -------
        s : `obj`:str:
            Multiline string describing this HDX Measurement object

        """

        timepoints = ", ".join([f"{t:.2f}" for t in self.timepoints])

        s = f"""
        HDX Measurement: {self.name}
        
        Number of peptides:     {self.Np}
        Number of residues:     {self.Nr} ({self.coverage.interval[0]} - {self.coverage.interval[1]})
        Number of timepoints:   {self.Nt}
        Timepoints:             {timepoints} seconds
        Coverage Percentage:    {self.coverage.percent_coverage:.2f}
        Average redundancy:     {self.coverage.redundancy:.2f}      
        Temperature:            {self.temperature} K
        pH:                     {self.pH}             
        """

        return textwrap.dedent(s.lstrip("\n"))

    def _repr_markdown_(self):
        s = str(self)
        s = s.replace("\n", "<br>")
        return s

    @property
    def name(self):
        """:obj:`str`: HDX Measurement name"""
        return self.metadata.get("name", self.state)

    @property
    def temperature(self):
        """:obj:`float`: Temperature of the H/D exchagne reaction (K)."""
        return self.metadata.get("temperature", None)

    @property
    def pH(self):
        """pH of the H/D exchange reaction."""
        return self.metadata.get("pH", None)

    @property
    def Np(self):
        """:obj:`int`: Number of peptides."""
        return self.coverage.Np

    @property
    def Nr(self):
        """:obj:`int`: Total number of residues spanned by the peptides."""
        return self.coverage.Nr

    @property
    def Nt(self):
        """:obj:`int`: Number of timepoints."""
        return len(self.timepoints)

    def __len__(self):
        import warnings

        warnings.warn("Use hdxm.Nt instead", DeprecationWarning)
        return len(self.timepoints)

    def __iter__(self):
        return self.peptides.__iter__()

    def __getitem__(self, item):
        return self.peptides.__getitem__(item)

    @property
    def rfu_residues(self):
        """:class:`~pandas.DataFrame`: Relative fractional uptake per residue. Shape Nr x Nt"""
        df = pd.concat([v.rfu_residues for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"

        return df

    @property
    def rfu_peptides(self):
        """:class:`~pandas.DataFrame`: Relative fractional uptake per peptide. Shape Np x Nt"""
        df = pd.concat([v.rfu_peptides for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"
        return df

    @property
    def d_exp(self):
        """:class:`~pandas.DataFrame`: D-uptake values (corrected). Shape Np x Nt"""
        df = pd.concat([v.d_exp for v in self], keys=self.timepoints, axis=1)
        df.columns.name = "exposure"
        return df

    def get_tensors(self, exchanges=False, dtype=None):
        """
        Returns a dictionary of tensor variables for fitting to Linderstrøm-Lang kinetics.

        Tensor variables are (shape):
        Temperature (1 x 1)
        X (Np x Nr)
        k_int (Nr x 1)
        timepoints (1 x Nt)
        d_exp (D) (Np x Nt)

        Parameters
        ----------
        exchanges : :obj:`bool`
            If True only returns tensor data describing residues which exchange (ie have peptides and are not prolines)

        Returns
        -------

        tensors : :obj:`dict`

        """

        if "k_int" not in self.coverage.protein:
            raise ValueError(
                "Unknown intrinsic rates of exchange, please supply pH and temperature parameters"
            )
        try:
            d_exp = self.d_exp
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

    def guess_deltaG(self, rates, correct_c_term=True):
        """
        Obtain ΔG initial guesses from apparent H/D exchange rates.
        Units of input rates are per second.

        Parameters
        ----------
        rates : :class:`~pandas.Series`
           Apparent exchange rate rates. Series index is protein residue number
        crop : :obj:`bool`
            If `True` the resulting :class:`~pandas.Series` is cropped to the residue interval covered by peptides.

        Returns
        -------
        dG : :class:`~pandas.Series`
            ΔG guess values

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
        file_path,
        include_version=True,
        include_metadata=True,
        fmt="csv",
        **kwargs,
    ):
        """
        Write the data in this HDX measurement to file.

        Parameters
        ----------
        file_path : :obj:`str`
            File path to create and write to.
        include_version : :obj:`bool`
            Set `True` to include PyHDX version and current time/date
        fmt: :obj: `str`
            Formatting to use, options are 'csv' or 'pprint'
        include_metadata : :obj:`bool`
            If `True`, the objects' metadata is included
        **kwargs : :obj:`dict`, optional
            Optional additional keyword arguments passed to `df.to_csv`
        Returns
        -------
        None

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
    """
    Class with subset of peptides corresponding to only one state and exposure

    Parameters
    ----------
    data : :class:`~pandas.DataFrame`
        Numpy structured array with input data

    """

    def __init__(self, data, **kwargs):
        assert len(np.unique(data["exposure"])) == 1, "Exposure entries are not unique"
        assert len(np.unique(data["state"])) == 1, "State entries are not unique"

        super(HDXTimepoint, self).__init__(data, **kwargs)

        self.state = self.data["state"][0]
        self.exposure = self.data["exposure"][0]

    @property
    def rfu_peptides(self):
        """:class:`~pandas.Series`: Relative fractional uptake per peptide"""
        return self.data["rfu"]

    @property
    def d_exp(self):
        """:class:`~pandas.Series`: Experimentally measured D-values (corrected)"""
        return self.data["uptake_corrected"]

    @property
    def name(self):
        """:obj:`str`: Name of this peptidemeasurement"""
        return self.state + "_" + str(self.exposure)

    @property
    def rfu_residues(self):
        """:class:`~pandas.Series`: Relative fractional uptake (RFU) per residue. Obtained by weighted averaging"""
        return self.weighted_average("rfu")

    def calc_rfu(self, residue_rfu):
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

    def weighted_average(self, field):
        """
        Calculate per-residue weighted average of values in data column

        Parameters
        ----------
        field : :obj:`str`
            Data field (column) to calculated weighted average of

        Returns
        -------


        """

        array = self.Z_norm.T.dot(self.data[field])
        series = pd.Series(array, index=self.index)
        return series


class CoverageSet(object):
    # todo perhaps this object should have X
    def __init__(self, hdxm_list):
        self.hdxm_list = hdxm_list

        # todo create Coverage object for the 3d case
        intervals = np.array(
            [hdxm_list.coverage.interval for hdxm_list in self.hdxm_list]
        )
        self.interval = (intervals[:, 0].min(), intervals[:, 1].max())
        self.r_number = np.arange(*self.interval)

        self.Ns = len(self.hdxm_list)
        self.Nr = len(self.r_number)
        self.Np = np.max([hdxm.Np for hdxm in self.hdxm_list])
        self.Nt = np.max([hdxm.Nt for hdxm in self.hdxm_list])

    @property
    def index(self):
        """pd index:"""
        return pd.RangeIndex(self.interval[0], self.interval[1], name="r_number")

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
    def s_r_mask(self):
        """mask of shape NsxNr with True entries covered by hdx measurements (exluding gaps)"""
        mask = np.zeros((self.Ns, self.Nr), dtype=bool)
        for i, hdxm in enumerate(self.hdxm_list):
            interval_sample = hdxm.coverage.interval
            i0 = interval_sample[0] - self.interval[0]
            i1 = interval_sample[1] - self.interval[0]

            mask[i, i0:i1] = True

        return mask

    def get_masks(self):
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
    Set of multiple :class:`~pyhdx.models.HDXMeasurement`

    Parameters
    ----------
    hdxm_list :  :obj:`list`
        or list of :class:`~pyhdx.models.HDXMeasurement`

    Attributes
    ----------
    timepoints : :class:`~numpy.ndarray`
        Ns x Nt array of zero-padded timepoints
    d_exp : :class:`~numpy.ndarray`
        Ns x Np x Nt array with zero-padded measured D-uptake values
    """

    def __init__(self, hdxm_list):
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

    def __getitem__(self, item):
        return self.hdxm_list.__getitem__(item)

    @property
    def Ns(self):
        return len(self.hdxm_list)

    @property
    def Nr(self):
        return self.coverage.Nr

    @property
    def Np(self):
        return np.max([hdxm.Np for hdxm in self.hdxm_list])

    @property
    def Nt(self):
        return np.max([hdxm.Nt for hdxm in self.hdxm_list])

    @property
    def temperature(self):
        return np.array([hdxm.temperature for hdxm in self.hdxm_list])

    @property
    def names(self):
        return [hdxm.name for hdxm in self.hdxm_list]

    @property
    def rfu_residues(self):
        rfu = pd.concat([hdxm.rfu_residues for hdxm in self],
                  keys=self.names, names=['state', 'exposure'], axis=1)
        columns = pd.MultiIndex.from_tuples(
            tuples=[(*tup, 'rfu') for tup in rfu.columns],
            names=['state', 'exposure', 'quantity']
        )
        
        rfu.columns = columns

        return rfu

    def guess_deltaG(self, rates_df):
        """
        Create dG guesses from rates

        Parameters
        ----------
        rates_df : :class:`~pandas.DataFrame`
            Pandas dataframe with k_obs estimates. Column names must correspond to HDX measurement names.

        Returns
        -------

        deltaG: :class:`~pandas.DataFrame`
            ΔG guess values

        """

        guesses = [
            hdxm.guess_deltaG(rates_df[name]) for hdxm, name in zip(self, self.names)
        ]
        deltaG = pd.concat(guesses, keys=self.names, axis=1)

        return deltaG

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

    def get_tensors(self, dtype=None):
        # todo create correct shapes as per table X for all
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
    def exchanges(self):
        values = np.concatenate(
            [hdxm.coverage["exchanges"].to_numpy() for hdxm in self.hdxm_list]
        )
        exchanges = np.zeros((self.Ns, self.Nr), dtype=bool)
        exchanges[self.masks["sr"]] = values

        return exchanges

    def to_file(
        self,
        file_path,
        include_version=True,
        include_metadata=True,
        fmt="csv",
        **kwargs,
    ):
        """
        Write the data in this HDX measurement set to file.

        Parameters
        ----------
        file_path : :obj:`str`
            File path to create and write to.
        include_version : :obj:`bool`
            Set ``True`` to include PyHDX version and current time/date
        fmt: :obj: `str`
            Formatting to use, options are 'csv' or 'pprint'
        include_metadata : :obj:`bool`
            If `True`, the objects' metadata is included
        **kwargs : :obj:`dict`, optional
            Optional additional keyword arguments passed to `df.to_csv`
        Returns
        -------
        None

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


def hdx_intersection(hdx_list, fields=None):
    """
    Finds the intersection between peptides in :class:`~pydhx.models.HDXMeasurement` and returns new objects such that
    all peptides (coverage, exposure) between the measurements are identical.

    Optionally intersections by custom fields can be made.

    Parameters
    ----------
    hdx_list : :obj:`list`
        Input list of :class:`~pyhdx.models.HDXMeasurement`
    fields : :obj:`list`
        By which fields to take the intersections. Default is ['_start', '_end', 'exposure']

    Returns
    -------
    hdx_out : :obj:`list`
        Output list of :class:`~pyhdx.models.HDXMeasurement`
    """

    fields = fields or ["_start", "_end", "exposure"]

    full_arrays = [data_obj.full_data for data_obj in hdx_list]
    selected = array_intersection(full_arrays, fields=fields)

    hdx_out = [
        HDXMeasurement(data, **data_obj.metadata)
        for data, data_obj in zip(selected, hdx_list)
    ]
    return hdx_out


def array_intersection(arrays_list, fields):
    """
    Find and return the intersecting entries in multiple arrays.

    Parameters
    ----------
    arrays_list : :obj:`iterable`
        Iterable of input structured arrays
    fields : :obj:`iterable`
        Iterable of fields to use to decide if entires are intersecting

    Returns
    -------
    selected : :obj:`iterable`
        Output iterable of arrays with only intersecting entries.

    """
    intersection = reduce(np.intersect1d, [fields_view(d, fields) for d in arrays_list])
    selected = [
        elem[np.isin(fields_view(elem, fields), intersection)] for elem in arrays_list
    ]

    return selected
