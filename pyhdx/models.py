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
from pyhdx.support import reduce_inter, make_view, fields_view


def protein_wrapper(func, *args, **kwargs):
    metadata = kwargs.pop('metadata', {})
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
                raise ValueError(f"Invalid index type {type(self.df.index)} for supplied DataFrame, must be integer index")

        if not self.df.index.is_unique:
            raise ValueError("Protein dataframe indices must be unique")

        new_index = pd.RangeIndex(start=self.df.index.min(), stop=self.df.index.max() + 1, name='r_number')
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

    def to_file(self, file_path, include_version=True, include_metadata=True, fmt='csv', **kwargs):
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
        dataframe_to_file(file_path, self.df, include_version=include_version, include_metadata=metadata, fmt=fmt, **kwargs)

    def set_k_int(self, temperature, pH):
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

        if 'sequence' not in self:
            raise ValueError('No sequence data available to calculate intrinsic exchange rates.')

        sequence = list(self['sequence'])  # Includes 'X' padding at cterm if cterm > last peptide
        k_int = k_int_from_sequence(sequence, temperature, pH)

        self.df['k_int'] = k_int
        return np.array(k_int)

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
    Main peptide input object. The input numpy structured array `data` must have the following entires for each peptide:

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
    data : :class:`~numpy.ndarray`
        Numpy recarray with peptide entries.
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

    def __init__(self, data, drop_first=1, ignore_prolines=True, d_percentage=100., sort=True, remove_nan=True):
        assert np.all(data['start'] < data['end']), 'All `start` entries must be smaller than their `end` entries'
        assert 0 <= d_percentage <= 100., 'Deuteration percentage must be between 0 and 100'
        d_percentage /= 100.

        self.data = data.copy()
        if remove_nan:
            self.data = self.data[~np.isnan(self.data['uptake'])]
        if sort:
            self.data = np.sort(self.data, order=['start', 'end', 'sequence', 'exposure', 'state'])

        # Make backup copies of unmodified start, end and sequence fields before taking prolines and n terminal residues into account
        if not np.any(np.isin(['_start', '_end', '_sequence'], self.data.dtype.names)):
            self.data = append_fields(self.data, ['_start'], [self.data['start'].copy()], usemask=False)
            self.data = append_fields(self.data, ['_end'], [self.data['end'].copy()], usemask=False)
            self.data = append_fields(self.data, ['_sequence'], [self.data['sequence'].copy()], usemask=False)

        # Convert sequence to upper case if not so already
        self.data['sequence'] = [s.upper() for s in self.data['sequence']]
        # Mark ignored prolines with lower case letters
        if ignore_prolines:
            self.data['sequence'] = [s.replace('P', 'p') for s in self.data['sequence']]

        # Find the total number of n terminal / c_terminal residues to remove
        # Todo: edge cases such as pure prolines or overlap between c terminal prolines and drop_first section (issue 32)
        n_term = np.array([len(seq) - len(seq[drop_first:].lstrip('p')) for seq in self.data['sequence']])
        c_term = np.array([len(seq) - len(seq.rstrip('p')) for seq in self.data['sequence']])

        # Mark removed n terminal residues with lower case x
        self.data['sequence'] = ['x'*nt + s[nt:] for nt, s in zip(n_term, self.data['sequence'])]
        self.data['start'] += n_term
        self.data['end'] -= c_term

        ex_residues = np.array([len(s) - s.count('x') - s.count('p') for s in self.data['sequence']]) * d_percentage
        if 'ex_residues' not in self.data.dtype.names:
            self.data = append_fields(self.data, ['ex_residues'], [ex_residues], usemask=False)

    def __len__(self):
        return len(self.data)

    def groupby_state(self, **kwargs):
        """
        Groups measurements in the dataset by state and returns them in a dictionary as a
        :class:`~pyhdx.models.HDXMeasurement`.

        Returns
        -------
        out : :obj:`dict`
            Dictionary where keys are state names and values are :class:`~pyhdx.models.HDXMeasurement`.
        **kwargs
            Additional keyword arguments to be passed to the :class:`~pyhdx.models.HDXMeasurement`.
        """

        warnings.warn("Likely to be removed in future versions, use `get_state` instead", PendingDeprecationWarning)

        states = np.unique(self.data['state'])
        return {state: HDXMeasurement(self.data[self.data['state'] == state], **kwargs) for state in states}

    def get_state(self, state):
        """
        Returns entries in the table with state 'state'

        Parameters
        ----------
        state : :obj:`str`


        Returns
        -------

        """
        return np.ascontiguousarray(self.data[self.data['state'] == state])

    @staticmethod
    def isin_by_idx(array, test_array):
        """
        Checks if entries in `array` are in `test_array`, by `start` and `end` field values.

        Parameters
        ----------
        array : :class:`~numpy.ndarray`
            Numpy input structured array
        test_array : :class:`~numpy.ndarray`
            Numpy structured array to test againts

        Returns
        -------
        isin : :obj:`ndarray`, :obj:`bool`
            Boolean array of the same shape as `array` where entries are `True` if they are in `test_array`

        """

        test = make_view(test_array, ['start', 'end'], dtype=np.int32)
        full = make_view(array, ['start', 'end'], dtype=np.int32)

        # https://stackoverflow.com/questions/54828039/how-to-match-pairs-of-values-contained-in-two-numpy-arrays/54828333
        isin = (full[:, None] == test).all(axis=2).any(axis=1)
        return isin

    def set_backexchange(self, back_exchange):
        """
        Sets the normalized percentage of uptake through a fixed backexchange value for all peptides.

        Parameters
        ----------
        back_exchange :  :obj:`float`
            Percentage of back exchange

        """

        back_exchange /= 100
        rfu = self.data['uptake'] / ((1-back_exchange)*self.data['ex_residues'])

        uptake_corrected = self.data['uptake'] / (1 - back_exchange)

        self.data = append_fields(self.data, ['rfu', 'uptake_corrected'], data=[rfu, uptake_corrected], usemask=False)

    def set_control(self, control_1, control_0=None):
        """
        Apply a control dataset to this object. A `scores` attribute is added to the object by normalizing its uptake
        value with respect to the control uptake value to 100%. Entries which are in the measurement and not in the
        control or vice versa are deleted.
        Optionally, ``control_zero`` can be specified which is a dataset whose uptake value will be used to zero
        the uptake.

        #todo insert math

        Parameters
        ----------
        control_1 : :obj:`tuple`
            tuple with (`state`, `exposure`) for peptides to use for normalization
        control_0 : :obj:`tuple`, optional
            tuple with (`state`, `exposure`) for peptides to use for zeroing uptake values

        """

        control_1 = self.get_data(*control_1)

        if control_0 is None:
            control_0 = np.copy(control_1)
            control_0['uptake'] = 0
        else:
            control_0 = self.get_data(*control_0)

        #todo this should probably go to the log (but atm there isnt any for running without GUI)
        assert control_1.size > 0, f"No peptides found with state '{control_1[0]}' and exposure '{control_1[1]}'"
        assert control_0.size > 0, f"No peptides found with state '{control_0[0]}' and exposure '{control_0[1]}'"

        b_100 = self.isin_by_idx(self.data, control_1)
        b_0 = self.isin_by_idx(self.data, control_0)
        data_selected = self.data[np.logical_and(b_100, b_0)]

        # Control peptides corresponding to those peptides in measurement
        c_1_selected = control_1[self.isin_by_idx(control_1, data_selected)]
        c_0_selected = control_0[self.isin_by_idx(control_0, data_selected)]

        control_1_final = np.sort(c_1_selected, order=['start', 'end', 'sequence', 'exposure', 'state'])
        control_0_final = np.sort(c_0_selected, order=['start', 'end', 'sequence', 'exposure', 'state'])

        # Sort both datasets by starting index and then by sequence to make sure they are both equal
        data_final = np.sort(data_selected, order=['start', 'end', 'sequence', 'exposure', 'state'])

        #Apply controls for each sequence  (#todo there must be a better way to do this)
        rfu = np.zeros(len(data_final), dtype=float)
        uptake_corrected = np.zeros(len(data_final), dtype=float)
        for c_1, c_0 in zip(control_1_final, control_0_final):
            bs = data_final['start'] == c_1['start']
            be = data_final['end'] == c_1['end']
            b_all = np.logical_and(bs, be)
            uptake = data_final[b_all]['uptake']
            rfu[b_all] = (uptake - c_0['uptake']) / (c_1['uptake'] - c_0['uptake'])

            uptake_corrected[b_all] = (uptake / c_1['uptake']) * data_final[b_all]['ex_residues']

        if 'rfu' in data_final.dtype.names:
            data_final['rfu'] = rfu
        else:
            data_final = append_fields(data_final, 'rfu', data=rfu, usemask=False)

        if 'uptake_corrected' in data_final.dtype.names:
            data_final['uptake_corrected'] = uptake_corrected
        else:
            data_final = append_fields(data_final, 'uptake_corrected', data=uptake_corrected, usemask=False)

        self.data = data_final

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
        output_data : :class:`~numpy.ndarray`
            Numpy structured array with selected peptides
        """

        output_data = self.data[np.logical_and(self.data['state'] == state, self.data['exposure'] == exposure)]
        return output_data

    @property
    def states(self):
        """:class:`~numpy.ndarray` Array with unique states"""
        return np.unique(self.data['state'])

    @property
    def exposures(self):
        """:class:`~numpy.ndarray` Array with unique exposures"""
        return np.unique(self.data['exposure'])


class Coverage(object):
    """
    Object describing layout and coverage of peptides and generating the corresponding matrices. Peptides should all
    belong to the same state and have the same exposure time.

    Parameters
    ----------
    data : :class:`~numpy.ndarray`
        Numpy structured array with input peptides
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

    def __init__(self, data, c_term=0, n_term=1, sequence=''):
        assert len(np.unique(data['exposure'])) == 1, 'Exposure entries are not unique'
        assert len(np.unique(data['state'])) == 1, 'State entries are not unique'

        self.data = np.sort(data, order=['start', 'end'])

        start = np.min(self.data['_start'])
        end = np.max(self.data['_end'])  # exclusive end interval

        if n_term:
            start = min(start, n_term)
        if sequence and not c_term:
            c_term = len(sequence) + n_term - 1
        if c_term:
            if c_term + 1 < end:
                raise ValueError("HDX data extends beyond c_term number, check 'sequence' or 'c_term'")
            end = c_term + 1  # c_term is inclusive, therefore plus one
        r_number = np.arange(start, end)  # r_number spanning the full protein range, not just the covered range
        # Full sequence
        _seq = np.full_like(r_number, fill_value='X', dtype='U')  # Full sequence
        # Sequence with lower case letters for no coverage due to n_terminal residues or prolines
        seq = np.full_like(r_number, fill_value='X', dtype='U')
        for d in self.data[::-1]:
            i, j = np.searchsorted(r_number, [d['_start'], d['_end']])
            _seq[i:j] = [s for s in d['_sequence']]
            seq[i:j] = [s for s in d['sequence']]

        if sequence:
            for r, s1, s2 in zip(r_number, sequence, _seq):
                if s2 != 'X' and s1 != s2:
                    raise ValueError(
                        f"Mismatch in supplied sequence and peptides sequence at residue {r}, expected '{s2}', got '{s1}'")
            if len(sequence) != len(_seq):
                raise ValueError("Invalid length of supplied sequence. Please check 'n_term' and 'c_term' parameters")
            _seq = list(sequence)

        #todo check if this is always correctly determined (n terminal residues usw)
        exchanges = [s.isupper() and (s != 'X') for s in seq]  # Boolean array True if residue exchanges, full length
        coverage = seq != 'X'  # Boolean array for coverage
        dic = {'r_number': r_number, 'sequence': _seq, 'coverage': coverage, 'exchanges': exchanges}

        # Inclusive, exclusive interval of peptides coverage across the whole protein
        self.interval = (np.min(self.data['start']), np.max(self.data['end']))
        self.protein = Protein(dic, index='r_number')

        # matrix dimensions N_peptides N_residues, dtype for TF compatibility
        _exchanges = self['exchanges']  # Array only on covered part
        self.X = np.zeros((len(self.data), self.interval[1] - self.interval[0]), dtype=int)
        self.Z = np.zeros_like(self.X, dtype=float)
        for row, entry in enumerate(self.data):
            i0, i1 = np.searchsorted(self.r_number, (entry['start'], entry['end']))
            self.X[row][i0:i1] = 1
            self.Z[row][i0:i1] = _exchanges[i0:i1]

        self.Z = self.Z / self.data['ex_residues'][:, np.newaxis]

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        pd_series = self.protein[item]
        return self.apply_interval(pd_series)

    def apply_interval(self, array_or_series):
        """Given a Numpy array or Pandas series with a length equal to the full protein, returns the section of the array equal to the covered
        region. Returned series length is equal to number of columns in the X matrix

        """

        if isinstance(array_or_series, np.ndarray):
            series = pd.Series(array_or_series, index=self.protein.df.index)
            assert len(array_or_series) == len(self.protein)
        else:
            series = array_or_series

        # - 1 because interval is inclusive, exclusive and .loc slices inclusive, inclusive
        covered_slice = series.loc[self.interval[0]:self.interval[1] - 1]

        return covered_slice

    @property
    def percent_coverage(self):
        """:obj:`float`: Percentage of residues covered by peptides"""
        return 100*np.mean(self.protein['coverage'])

    @property
    def redundancy(self):
        """:obj:`float`: Average redundancy of peptides in regions with at least 1 peptide"""
        x_coverage = self.X[:, self['coverage']]
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
        """:class:`~numpy.ndarray`: Array of residue numbers corresponding to the part of the protein covered by peptides"""
        #todo perhaps obtain through apply_interval
        return np.arange(*self.interval)

    @property
    def index(self):
        """:class:`~pandas.RangeIndex` """
        return pd.RangeIndex(self.interval[0], self.interval[1], name='r_number')
        #return pd.Index(self.r_number, name='r_number')

    @property
    def block_length(self):
        """:class:`~numpy.ndarary`: Lengths of unique blocks of residues in the peptides map,
            along the `r_number` axis"""

        # indices are start and stop values of blocks
        indices = np.sort(np.concatenate([self.data['start'], self.data['end']]))
        #indices of insertion into r_number vector gives us blocks with taking prolines into account.
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
        return self.Z / np.sum(self.Z, axis=0)[np.newaxis, :]

    def get_sections(self, gap_size=-1):
        """get the intervals of sections of coverage
        intervals are inclusive, exclusive

            gap_size : :obj:`int`
        Gaps of this size between adjacent peptides is not considered to overlap. A value of -1 means that peptides
        with exactly zero overlap are separated. With gap_size=0 peptides with exactly zero overlap are not separated,
        and larger values tolerate larger gap sizes.
        """
        intervals = [(s, e) for s, e in zip(self.data['start'], self.data['end'])]
        sections = reduce_inter(intervals, gap_size=gap_size)

        return sections

    def __eq__(self, other):
        """Coverage objects are considered equal if both objects fully match between their start, end and sequence fields"""
        assert isinstance(other, Coverage), "Other must be an instance of Coverage"
        return len(self.data) == len(other.data) and np.all(self.data['start'] == other.data['start']) and \
               np.all(self.data['end'] == other.data['end']) and np.all(self.data['sequence'] == other.data['sequence'])


class HDXMeasurement(object):
    """
    Main HDX data object. This object has peptide data of a single state but with multiple timepoints.

    Timepoint data is split into :class:`~pyhdx.models.PeptideMeasurements` objects for each timepoint
    Supplied data is made 'uniform' such that all timepoints have the same peptides

    Parameters
    ----------
    data : :class:`~numpy.ndarray` or :obj:`list`
        Numpy structured array with peptide entries corresponding to a single state,
        or list of :class:`~pyhdx.models.PeptideMeasurements`
    **metadata
        Dictionary of optional metadata. By default, holds the `temperature` and `pH` parameters.


    Attributes
    ----------
    state : :obj:`str`
        State of the HDX measurement
    timepoints : :class:`~numpy.ndarray`
        Array with exposure times (sorted)
    peptides : :obj:`list`
        List of :class:`~pyhdx.models.PeptideMeasurements`, one list element per timepoint.
    cov : :class:`~pyhdx.models.Coverage`
        Coverage object describing peptide layout. When this `uniform` is `False`, this attribute is `None`

    """
    def __init__(self, data, **metadata):
        self.metadata = metadata
        assert len(np.unique(data['state'])) == 1
        self.state = str(data['state'][0])
        self.timepoints = np.sort(np.unique(data['exposure']))

        data_list = [(data[data['exposure'] == exposure]) for exposure in self.timepoints]
        sets = [{tuple(elem) for elem in fields_view(d, ['_start', '_end'])} for d in data_list]
        intersection = set.intersection(*sets)
        dtype = [('_start', data['_start'].dtype), ('_end', data['_end'].dtype)]
        intersection_array = np.array([tup for tup in intersection], dtype=dtype)

        # Select entries in data array which are in the intersection between all timepoints
        selected = [elem[np.isin(fields_view(elem, ['_start', '_end']), intersection_array)] for elem in data_list]
        cov_kwargs = {kwarg: metadata.get(kwarg, default) for kwarg, default in zip(['c_term', 'n_term', 'sequence'], [0, 1, ''])}
        self.peptides = [PeptideMeasurements(elem, **cov_kwargs) for elem in selected]

        # Create coverage object from the first time point (as all are now equal)
        self.coverage = Coverage(selected[0], **cov_kwargs)

        if self.temperature and self.pH:
            self.coverage.protein.set_k_int(self.temperature, self.pH)

    def __str__(self):
        """

        Returns
        -------
        s : `obj`:str:
            Multiline string describing this HDX Measurement object

        """

        timepoints = ', '.join([f'{t:.2f}' for t in self.timepoints])

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

        return textwrap.dedent(s)


    @property
    def name(self):
        return self.metadata.get('name', self.state)

    @property
    def temperature(self):
        return self.metadata.get('temperature', None)

    @property
    def pH(self):
        return self.metadata.get('pH', None)

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

    @property
    def full_data(self):
        """returns the full dataset of all timepoints"""
        full_data = np.concatenate([pm.data for pm in self])
        return full_data

    def __len__(self):
        return len(self.timepoints)

    def __iter__(self):
        return self.peptides.__iter__()

    def __getitem__(self, item):
        return self.peptides.__getitem__(item)

    @property
    def rfu_residues(self):
        """Relative fractional uptake per residue. Shape Nr x Nt"""
        return np.stack([v.rfu_residues for v in self]).T

    @property
    def rfu_peptides(self):
        return np.stack([v.rfu_peptides for v in self])

    @property
    def uptake_corrected(self):
        """matrix shape  N_t, N_p""" #(should be np nt)
        #todo refactor to D to match manuscript
        #todo deprecate
        uptake_corrected = np.stack([v.uptake_corrected for v in self])
        return uptake_corrected

    @property
    def d_exp(self):
        """np.ndarray (shape Np x Nt)
            Experimentally measured D-uptake values, corrected for back-exchange """
        return self.uptake_corrected.T

    def get_tensors(self, exchanges=False):
        """
        Returns a dictionary of tensor variables for fitting to Linderstrøm-Lang kinetics.

        Tensor variables are (shape):
        Temperature (1 x 1)
        X (Np x Nr)
        k_int (Nr x 1)
        timepoints (1 x Nt)
        uptake (D) (Np x Nt)

        Parameters
        ----------
        exchanges : :obj:`bool`
            if True only returns tensor data describing residues which exchange (ie have peptides and are not prolines)

        Returns
        -------

        tensors : :obj:`dict`

        """

        if 'k_int' not in self.coverage.protein:
            raise ValueError("Unknown intrinsic rates of exchange, please supply pH and temperature parameters")
        try:
            upt = self.uptake_corrected
        except ValueError:
            raise ValueError("HDX data is not corrected for back exchange.")

        if exchanges:
            #this could be a method on coverage object similar to apply_interval; select exchanging
            bools = self.coverage['exchanges'].to_numpy()
        else:
            bools = np.ones(self.Nr, dtype=bool)

        dtype = pyhdx.fitting_torch.TORCH_DTYPE
        device = pyhdx.fitting_torch.TORCH_DEVICE

        tensors = {
            'temperature': torch.tensor([self.temperature], dtype=dtype, device=device).unsqueeze(-1),
            'X': torch.tensor(self.coverage.X[:, bools], dtype=dtype, device=device),
            'k_int': torch.tensor(self.coverage['k_int'].to_numpy()[bools], dtype=dtype, device=device).unsqueeze(-1),
            'timepoints': torch.tensor(self.timepoints, dtype=dtype, device=device).unsqueeze(0),
            'uptake': torch.tensor(self.uptake_corrected.T, dtype=dtype, device=device)}

        return tensors

    def guess_deltaG(self, rates, crop=True):
        """

        Parameters
        ----------
        rates : :class:`~pandas.Series`
            pandas series of estimated hdx exchangs rates. Index is protein residue number
        return_type

        Returns
        -------

        """
        if 'k_int' not in self.coverage.protein:
            raise ValueError("Unknown intrinsic rates of exchange, please supply pH and temperature parameters")
        if not isinstance(rates, pd.Series):
            raise TypeError("Rates input type is pandas.Series")

        p_guess = (self.coverage.protein['k_int'] / rates) - 1

        p_guess.clip(0., None, inplace=True)  # Some initial guesses might have negative PF values
        with np.errstate(divide='ignore'):
            deltaG = np.log(p_guess) * constants.R * self.temperature

        # https://stackoverflow.com/questions/9537543/replace-nans-in-numpy-array-with-closest-non-nan-value
        bools = ~np.isfinite(deltaG)
        deltaG[bools] = np.interp(np.flatnonzero(bools), np.flatnonzero(~bools), deltaG[~bools])

        if crop:
            return self.coverage.apply_interval(deltaG)
        else:
            return deltaG

    def to_file(self, file_path, include_version=True, include_metadata=True, fmt='csv', **kwargs):
        """
        Write the data in this HDX measurement to file.


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
        #todo save full_data as pandas dataframe

        # requires testing dont think this works as intended
        # should use self.metadata if include_metadata is the bool True otherwise if its a dict use that
        metadata = self.metadata if include_metadata else include_metadata
        df = pd.DataFrame(self.full_data)
        df.index.name = 'peptide_index'
        df.index += 1
        dataframe_to_file(file_path, df, include_version=include_version, include_metadata=metadata, fmt=fmt, **kwargs)


class PeptideMeasurements(Coverage):
    """
    Class with subset of peptides corresponding to only one state and exposure

    Parameters
    ----------
    data : :class:`~numpy.ndarray`
        Numpy structured array with input data
    scores : :class:`~numpy.ndarray`
        Array with D/H uptake scores, typically in percentages or absolute uptake numbers.

    Attributes
    ----------
    start : :obj:`int`
        First peptide starts at this residue number (starting from 1)
    stop : :obj:`int`
        Last peptide ends at this residue number (incusive)
    prot_len : :obj:`int`
        Total number of residues in this set of peptides, not taking regions of no coverage into account.
    exposure : :obj:`float`
        Exposure time of this set of peptides (minutes)
    state : :obj:`string`
        State describing the experiment

    bigX
    X

    properties:
    big_x_norm
    x_norm

    scores nnls
    scores lsq

    """

    def __init__(self, data, **kwargs):
        assert len(np.unique(data['exposure'])) == 1, 'Exposure entries are not unique'
        assert len(np.unique(data['state'])) == 1, 'State entries are not unique'

        super(PeptideMeasurements, self).__init__(data, **kwargs)

        self.state = self.data['state'][0]
        self.exposure = self.data['exposure'][0]

    @property
    def rfu_peptides(self):
        return self.data['rfu']

    @property
    def uptake_corrected(self):
        return self.data['uptake_corrected']

    @property
    def name(self):
        return self.state + '_' + str(self.exposure)

    @property
    def rfu_residues(self):
        """Weighted averaged relative fractional uptake"""
        return self.Z_norm.T.dot(self.rfu_peptides)

    def calc_rfu(self, residue_rfu):
        """
        Calculates RFU per peptide given an array of individual residue scores

        Parameters
        ----------
        residue_rfu : :class:`~numpy.ndarray`
            Array of rfu per residue of length `prot_len`

        Returns
        -------

        scores : :class:`~numpy.ndarray`
            Array of rfu per peptide
        """

        scores = self.Z.dot(residue_rfu)
        return scores

    def weighted_average(self, field):
        """Calculate per-residue weighted average of values in data column given by 'field'"""

        return self.Z_norm.T.dot(self.data[field])


class CoverageSet(object):
    #todo perhaps this object should have X
    def __init__(self, hdxm_list):
        self.hdxm_list = hdxm_list

        #todo create Coverage object for the 3d case
        intervals = np.array([hdxm_list.coverage.interval for hdxm_list in self.hdxm_list])
        self.interval = (intervals[:, 0].min(), intervals[:, 1].max())
        self.r_number = np.arange(*self.interval)

        self.Ns = len(self.hdxm_list)
        self.Nr = len(self.r_number)
        self.Np = np.max([hdxm.Np for hdxm in self.hdxm_list])
        self.Nt = np.max([hdxm.Nt for hdxm in self.hdxm_list])

    @property
    def index(self):
        """pd index: """
        return pd.RangeIndex(self.interval[0], self.interval[1], name='r_number')

    def apply_interval(self, array_or_series):
        """Given a Numpy array or Pandas series with a length equal to the full protein, returns the section of the array equal to the covered
        region. Returned series length is equal to number of columns in the X matrix

        """
        #todo testing and 2d array support
        if isinstance(array_or_series, np.ndarray):
            series = pd.Series(array_or_series, index=self.index)
            assert len(array_or_series) == len(self.index)
        else:
            series = array_or_series

        # - 1 because interval is inclusive, exclusive and .loc slices inclusive, inclusive
        covered_slice = series.loc[self.interval[0]:self.interval[1] - 1]

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
            st_mask[i, -hdxm.Nt:] = True
            spr_mask[i, 0: hdxm.Np, i0:i1] = True
            spt_mask[i, 0: hdxm.Np, -hdxm.Nt:] = True

        mask_dict = {'sr': sr_mask, 'st': st_mask, 'spr': spr_mask, 'spt': spt_mask}

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
        self.timepoints[self.masks['st']] = timepoints_values

        d_values = np.concatenate([hdxm.d_exp.flatten() for hdxm in self.hdxm_list])
        self.d_exp = np.zeros((self.Ns, self.Np, self.Nt))
        self.d_exp[self.masks['spt']] = d_values

        # Index array of of shape Ns x y where indices apply to deltaG return aligned residues for
        self.aligned_indices = None
        self.aligned_dataframes = None

    def __iter__(self):
        return self.hdxm_list.__iter__()

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

    def guess_deltaG(self, rates_list):
        """
        create deltaG guesses from rates

        Parameters
        ----------
        rates_list : :obj:`iterable`
            list of pandas series with k_obs esimates

        Returns
        -------

        deltaG_array: :class:`~numpy.ndarray`
            deltaG guesses Ns x Nr shape

        """
        assert len(rates_list) == self.Ns, "Number of elements in 'rates_list' should be equal to number of samples"

        guesses = [hdxm.guess_deltaG(rates, crop=True).to_numpy() for rates, hdxm in zip(rates_list, self.hdxm_list)]
        flat = np.concatenate(guesses)

        deltaG_array = np.full((self.Ns, self.Nr), fill_value=np.nan)
        deltaG_array[self.coverage.s_r_mask] = flat  # todo get this mask from dict?

        for row in deltaG_array:
            # https://stackoverflow.com/questions/9537543/replace-nans-in-numpy-array-with-closest-non-nan-value
            bools = ~np.isfinite(row)
            row[bools] = np.interp(np.flatnonzero(bools), np.flatnonzero(~bools), row[~bools])

        return deltaG_array

    def add_alignment(self, alignment, first_r_numbers=None):
        """

        :param alignment: list
        :param first_r_numbers:
            default is [1, 1, ...] but specifiy here if alignments do not all start at residue 1
        :return:
        """
        dfs = [hdxm.coverage.protein.df for hdxm in self.hdxm_list]
        self.aligned_dataframes = align_dataframes(dfs, alignment, first_r_numbers)

        df = self.aligned_dataframes['r_number']

        # Crop residue numbers to interval range
        df = df[((self.coverage.interval[0] <= df) & (df < self.coverage.interval[1])).all(axis=1)]
        df = df - self.coverage.interval[0]  # First residue in interval selected by index 0
        df.dropna(how='any', inplace=True)  # Remove non-aligned residues

        self.aligned_indices = df.to_numpy(dtype=int).T

    def get_tensors(self):
        #todo create correct shapes as per table X for all
        temperature = np.array([kf.temperature for kf in self.hdxm_list])

        X_values = np.concatenate([hdxm.coverage.X.flatten() for hdxm in self.hdxm_list])
        X = np.zeros((self.Ns, self.Np, self.Nr))
        X[self.masks['spr']] = X_values

        k_int_values = np.concatenate([hdxm.coverage['k_int'].to_numpy() for hdxm in self.hdxm_list])
        k_int = np.zeros((self.Ns, self.Nr))
        k_int[self.masks['sr']] = k_int_values

        dtype = pyhdx.fitting_torch.TORCH_DTYPE
        device = pyhdx.fitting_torch.TORCH_DEVICE

        tensors = {
            'temperature': torch.tensor(temperature, dtype=dtype, device=device).reshape(self.Ns, 1, 1),
            'X': torch.tensor(X, dtype=dtype, device=device),
            'k_int': torch.tensor(k_int, dtype=dtype, device=device).reshape(self.Ns, self.Nr, 1),
            'timepoints': torch.tensor(self.timepoints, dtype=dtype, device=device).reshape(self.Ns, 1, self.Nt),
            'uptake': torch.tensor(self.d_exp, dtype=dtype, device=device)  #todo this is called uptake_corrected/D/uptake
        }

        return tensors

    @property
    def exchanges(self):
        values = np.concatenate([hdxm.coverage['exchanges'].to_numpy() for hdxm in self.hdxm_list])
        exchanges = np.zeros((self.Ns, self.Nr), dtype=bool)
        exchanges[self.masks['sr']] = values

        return exchanges

    def to_file(self, file_path, include_version=True, include_metadata=True, fmt='csv', **kwargs):
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
            metadata[hdxm.name] = hdxm.metadata if include_metadata else include_metadata
            df = pd.DataFrame(hdxm.full_data)
            df.index.name = 'peptide_index'
            df.index += 1
            dfs.append(df)

        full_df = pd.concat(dfs, axis=1, keys=self.names)
        dataframe_to_file(file_path, full_df, include_version=include_version, include_metadata=metadata, fmt=fmt, **kwargs)


#https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)
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

    fields = fields or ['_start', '_end', 'exposure']

    full_arrays = [data_obj.full_data for data_obj in hdx_list]
    selected = array_intersection(full_arrays, fields=fields)

    hdx_out = [HDXMeasurement(data, **data_obj.metadata) for data, data_obj in zip(selected, hdx_list)]
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
    selected = [elem[np.isin(fields_view(elem, fields), intersection)] for elem in arrays_list]

    return selected
