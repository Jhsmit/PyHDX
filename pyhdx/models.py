import numpy as np
from numpy.lib.recfunctions import append_fields
import itertools
import scipy
from datetime import datetime
import pandas as pd
from io import StringIO
from functools import reduce
from operator import add
from hdxrate import k_int_from_sequence
from pyhdx.support import reduce_inter, make_view, fields_view
from pyhdx.fileIO import fmt_export
import pyhdx


class Protein(object):
    """Object describing a protein

    Parameters
    ----------
    data : :class:`~np.ndarray` or ?
        data object to initiate the protein object from
    index: :obj:`str`
        Name of the column with the residue number (index column)

    **metadata
        Dictionary of optional metadata.


    """

    def __init__(self, data, index, **metadata):
        self.metadata = metadata
        self.df = pd.DataFrame(data)
        if self.df.index.name is None:
            self.df.set_index(index, inplace=True)

        self.df.sort_values(index, inplace=True)

    def __str__(self):
        s = self.df.__str__()
        full_s = "Protein <name>\n" + s
        return full_s

    def __len__(self):
        return len(self.df)

    def join(self, other, on=None, how='left', lsuffix='', rsuffix='', sort=False):
        """
        Metadata is merged (overlapping values are taken from other)

        Parameters
        ----------
        other
        on
        how
        lsuffix
        rsuffix
        sort

        Returns
        -------

        """
        df_out = self.df.join(other.df, on=on, how=how, lsuffix=lsuffix, rsuffix=rsuffix, sort=sort)
        return self._make_protein(df_out, other)

    def append(self, other, ignore_index=False, verify_integrity=False, sort=None):
        df_out = self.df.append(other.df, ignore_index=ignore_index, verify_integrity=verify_integrity, sort=sort)
        return self._make_protein(df_out, other)

    def concat(self, other):
        df_out = pd.concat([self.df, other.df], axis=1)
        return self._make_protein(df_out, other)

    def merge(self, other, ignore_index=False, verify_integrity=False, sort=None):
        df_out = self.df.merge(other.df, how='left', validate='m:m')
        return self._make_protein(df_out, other)

    def add_column(self, index, data, name):
        pass

    def _make_protein(self, df_out, other):
        """Make a new :class:`~pyhdx.models.Protein` object and combine metadata with other metadata"""
        metadata = {**self.metadata, **other.metadata}
        protein_out = Protein(df_out, index=df_out.index.name, **metadata)
        return protein_out

    def to_records(self, index=True, column_dtypes=None, index_dtypes=None):
        return self.df.to_records(index=index, column_dtypes=column_dtypes, index_dtypes=index_dtypes)

    def to_stringio(self, io=None, include_version=True, include_metadata=True):
        """
        Write Protein data to :class:`~io.StringIO`

        Parameters
        ----------
        io : :class:`~io.StringIO`, optional
            StringIO to write to. If `None` a new StringIO object is created.
        include_version : :obj:`bool`
            Set `True` to include PyHDX version and current time/date
        include_metadata
            Not Implemented

        Returns
        -------
        io : :class:`~io.StringIO`
        """
        #todo add metadata

        io = io or StringIO()

        if include_version:
            io.write('# ' + pyhdx.VERSION_STRING + ' \n')
            now = datetime.now()
            io.write(f'# {now.strftime("%Y/%m/%d %H:%M:%S")} ({int(now.timestamp())}) \n')

        records = self.to_records()
        fmt, hdr = fmt_export(records)
        np.savetxt(io, records, fmt=fmt, header=hdr)

        io.seek(0)
        return io

    def to_file(self, file_path, include_version=True, include_metadata=True):
        """
        Write Protein data to file.


        Parameters
        ----------
        file_path : :obj:`str`
            File path to create and write to.
        include_version : :obj`bool`
            Set `True` to include PyHDX version and current time/date
        include_metadata
            Not Implemented

        Returns
        -------

        None

        """
        io = self.to_stringio(include_version=include_version, include_metadata=include_metadata)
        with open(file_path, 'w') as f:
            print(io.getvalue(), file=f)

    def to_dict(self, orient='list', into=dict):
        """
        Panda's `to_dict`, with different default options

        Parameters
        ----------
        orient
        into

        Returns
        -------

        """
        return self.df.to_dict(orient=orient, into=into)

    def set_k_int(self, temperature, pH):
        """
        Calculates the intrinsic rate of the sequence. Values of no coverage or prolines are assigned a value of -1
        The rates run are for the first residue (1) up to the last residue that is covered by peptides

        When the previous residue is unknown the current residue is also assigned a value of -1.g

        Parameters
        ----------
        temperature: : :obj:`float`
            Temperature of the labelling reaction (Kelvin)
        pH : :obj:`float`
            pH of the labelling reaction

        Returns
        -------

        k_int : ~class:`~numpy.ndarray`
            Array of intrisic exchange rates

        """

        if 'sequence' not in self:
            raise ValueError('No sequence data available to calculate intrinsic exchange rates.')

        sequence = list(self['sequence'])  # Includes 'X' padding at cterm if cterm > last peptide
        k_int = k_int_from_sequence(sequence, temperature, pH) * 60  # convert to per minute from per second

        self.df['k_int'] = k_int
        return np.array(k_int)

    @property
    def c_term(self):
        return self.df.index.max()

    @property
    def index(self):
        return self.df.index

    def copy(self):
        df = self.df.copy()
        return Protein(df, index=df.index.name, **self.metadata.copy())  #todo

    def __getitem__(self, item):
        return self.df.__getitem__(item)

    def __contains__(self, item):
        return self.df.__contains__(item)

    def __sub__(self, other):
        assert isinstance(other, Protein)
        df_out = self.df.subtract(other.df)
        return self._make_protein(df_out, other)

    def __add__(self, other):
        assert isinstance(other, Protein)
        df_out = self.df.add(other.df)
        return self._make_protein(df_out, other)

    def __truediv__(self, other):
        assert isinstance(other, Protein)
        df_out = self.df.truediv(other.df)
        return self._make_protein(df_out, other)

    def __floordiv__(self, other):
        assert isinstance(other, Protein)
        df_out = self.df.floordiv(other.df)
        return self._make_protein(df_out, other)

    def __mul__(self, other):
        assert isinstance(other, Protein)
        df_out = self.df.mul(other.df)
        return self._make_protein(df_out, other)


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
    data : ~:class:`np.ndarray`
        Numpy recarray with peptide entries.
    drop_first : :obj:`int`
        Number of N-terminal amino acids to ignore. Default is 1.
    d_percentage : :obj:`float`
        Percentage of deuterium in the labelling solution.
    ignore_prolines: :obj:`bool`
        Boolean to toggle ignoring of proline residues. When True these residues are treated as if they're not present
        in the protein.
    sort: :obj:`bool`
        Set to ``True`` to sort the input. Sort order is 'start', 'end', 'sequence', 'exposure', 'state'.
    remove_nan: :obj`bool`
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
        :class:`~pyhdx.pyhdx.KineticSeries`.

        Returns
        -------
        out : :obj:`dict`
            Dictionary where keys are state names and values are :class:`~pyhdx.pyhdx.KineticSeries`.
        **kwargs
            Additional keyword arguments to be passed to the :class:`~pyhdx.pyhdx.KineticSeries`.
        """

        states = np.unique(self.data['state'])
        return {state: KineticsSeries(self.data[self.data['state'] == state], **kwargs) for state in states}

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
        isin: ndarray, bool
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
        back_exchange :  `obj`:float:
            Percentage of back exchange

        """

        back_exchange /= 100
        scores = 100*self.data['uptake'] / ((1-back_exchange)*self.data['ex_residues'])

        uptake_corrected = self.data['uptake'] / (1 - back_exchange)

        self.data = append_fields(self.data, ['scores', 'uptake_corrected'], data=[scores, uptake_corrected], usemask=False)

    def set_control(self, control_100, control_0=None):
        """
        Apply a control dataset to this object. A `scores` attribute is added to the object by normalizing its uptake
        value with respect to the control uptake value to 100%. Entries which are in the measurement and not in the
        control or vice versa are deleted.
        Optionally, ``control_zero`` can be specified which is a dataset whose uptake value will be used to zero
        the uptake.

        #todo insert math

        Parameters
        ----------
        control_100 : tuple
            tuple with (`state`, `exposure`) for peptides to use for normalization to 100%
            Numpy structured array with control peptides to use for normalization to 100%
        control_0 : tuple, optional
            tuple with (`state`, `exposure`) for peptides to use for zeroing uptake values to 100%

        """

        control_100 = self.get_data(*control_100)

        if control_0 is None:
            control_0 = np.copy(control_100)
            control_0['uptake'] = 0
        else:
            control_0 = self.get_data(*control_0)

        b_100 = self.isin_by_idx(self.data, control_100)
        b_0 = self.isin_by_idx(self.data, control_0)
        data_selected = self.data[np.logical_and(b_100, b_0)]

        # Control peptides corresponding to those peptides in measurement
        c_100_selected = control_100[self.isin_by_idx(control_100, data_selected)]
        c_0_selected = control_0[self.isin_by_idx(control_0, data_selected)]

        control_100_final = np.sort(c_100_selected, order=['start', 'end', 'sequence', 'exposure', 'state'])
        control_0_final = np.sort(c_0_selected, order=['start', 'end', 'sequence', 'exposure', 'state'])

        # Sort both datasets by starting index and then by sequence to make sure they are both equal
        data_final = np.sort(data_selected, order=['start', 'end', 'sequence', 'exposure', 'state'])

        #Apply controls for each sequence  (#todo there must be a better way to do this)
        scores = np.zeros(len(data_final), dtype=float)
        uptake_corrected = np.zeros(len(data_final), dtype=float)
        for c_100, c_0 in zip(control_100_final, control_0_final):
            bs = data_final['start'] == c_100['start']
            be = data_final['end'] == c_100['end']
            b_all = np.logical_and(bs, be)
            uptake = data_final[b_all]['uptake']
            scores[b_all] = 100 * (uptake - c_0['uptake']) / (c_100['uptake'] - c_0['uptake'])

            uptake_corrected[b_all] = (uptake / c_100['uptake']) * data_final[b_all]['ex_residues']

        if 'scores' in data_final.dtype.names:
            data_final['scores'] = scores
        else:
            data_final = append_fields(data_final, 'scores', data=scores, usemask=False)

        if 'uptake_corrected' in data_final.dtype.names:
            data_final['uptake_corrected'] = uptake_corrected
        else:
            data_final = append_fields(data_final, 'uptake_corrected', data=uptake_corrected, usemask=False)

        self.data = data_final

    def return_by_name(self, control_state, control_exposure):
        #todo return dictionary of kinetic series instead
        print('deprecate this')  #currently used by GUI
        #raise DeprecationWarning

        """

        Finds all peptides in the dataset which match the control peptides and the peptides are grouped by their state
        and exposure and returned in a dictionary.

        Parameters
        ----------
        control_state : :obj:`str`
            Name of the control state
        control_exposure : :obj:`float`
            Exposure time of the control

        Returns
        -------
        out : :obj:`dict`
            Dictionary of :class:`~pyhdx.pyhdx.PeptideMeasurement` objects

        """
        bc1 = self.data['state'] == control_state
        bc2 = self.data['exposure'] == control_exposure

        control = self.data[np.logical_and(bc1, bc2)]
        st = np.unique(self.data['state'])
        exp = np.unique(self.data['exposure'])

        out = {}
        # Iterative over all permutations of state and exposure time to find all entries
        for s, e in itertools.product(st, exp):
            b1 = self.data['state'] == s
            b2 = self.data['exposure'] == e
            bf = np.logical_and(b1, b2)
            name = s + '_' + str(round(e, 3))

            d = self.data[bf]
            b_data = np.isin(d['sequence'], control['sequence'])  # find the sequences in the measurement that are in control
            d_selected = d[b_data]
            b_control = np.isin(control['sequence'], d_selected['sequence']) # find the control entries corresponding to these sequences
            control_selected = control[b_control]

            #sort both datasets by starting index then by sequence
            data_final = np.sort(d_selected, order=['start', 'sequence'])
            control_final = np.sort(control_selected, order=['start', 'sequence'])

            assert np.all(data_final['sequence'] == control_final['sequence'])
            assert np.all(data_final['start'] == control_final['start'])
            assert np.all(data_final['end'] == control_final['end'])
            score = 100 * data_final['uptake'] / control_final['uptake']

            if len(score) > 0:
                out[name] = PeptideMeasurements(data_final)

        return out

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
        """:~classs:np.ndarray: Array with unique states"""
        return np.unique(self.data['state'])

    @property
    def exposures(self):
        """:~classs:np.ndarray: Array with unique exposures"""
        return np.unique(self.data['exposure'])


class Coverage(object):
    """
    Object describing layout and coverage of peptides and generating the corresponding matrices. Peptides should all
    belong to the same state and have the same exposure time.

    Parameters
    ----------
    data : ~class:`~numpy.ndarray`
        Numpy structured array with input peptides

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

    def __init__(self, data, c_term=None):
        assert len(np.unique(data['exposure'])) == 1, 'Exposure entries are not unique'
        assert len(np.unique(data['state'])) == 1, 'State entries are not unique'

        self.data = np.sort(data, order=['start', 'end'])

        start = min(np.min(self.data['_start']), 1)
        end = np.max(self.data['_end'])
        if c_term:
            end = max(end, c_term + 1)  # c_term is inclusive, therefore plus one
        r_number = np.arange(start, end)  # r_number spanning the full protein range, not just the covered range

        # Full sequence
        _seq = np.full_like(r_number, fill_value='X', dtype='U')  # Full sequence
        # Sequence with lower case letters for no coverage due to n_terminal residues or prolines
        seq = np.full_like(r_number, fill_value='X', dtype='U')
        for d in self.data[::-1]:
            i, j = np.searchsorted(r_number, [d['_start'], d['_end']])
            _seq[i:j] = [s for s in d['_sequence']]
            seq[i:j] = [s for s in d['sequence']]

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
        series = self.protein[item]
        #return self.apply_interval(series.to_numpy())
        return self.apply_interval(series)

    def apply_interval(self, array_or_series):
        """Given an array or series with a length equal to the full protein, returns the section of the array equal to the covered
        region. Returned series length is equal to number of colunms in the X matrix

        """

        assert len(array_or_series) == len(self.protein)
        if isinstance(array_or_series, np.ndarray):
            series = pd.Series(array_or_series, index=self.protein.df.index)
        else:
            series = array_or_series

        # - 1 because interval is inclusive, exclusive and .loc slices inclusive, inclusive
        covered_slice = series.loc[self.interval[0]:self.interval[1] - 1]

        return covered_slice

    @property
    def r_number(self):
        """:class:`~np.ndarray: Array of residue numbers corresponding to the part of the protein covered by peptides"""
        #todo perhaps obtain through apply_interval
        return np.arange(*self.interval)

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
        """:class:`~np.ndarray`: `X` coefficient matrix normalized column wise."""
        return self.X / np.sum(self.X, axis=0)[np.newaxis, :]

    @property
    def Z_norm(self):
        """:class:`~np.ndarray`: `Z` coefficient matrix normalized column wise."""
        return self.Z / np.sum(self.Z, axis=0)[np.newaxis, :]

    def get_sections(self, gap_size=-1):
        """get the intervals of sections of coverage
        intervals are inclusive, exclusive

            gap_size: :obj:`int`
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


class KineticsSeries(object):
    """
    A series of :class:`~pyhdx.pyhdx.PeptideMeasurements` which correspond to the same state but with different exposures.

    Parameters
    ----------
    data : :class:`~numpy.ndarray` or :obj:`list`
        Numpy structured array with peptide entries corresponding to a single state,
        or list of :class:`~pyhdx.pyhdx.PeptideMeasurements`
    make_uniform : :obj:`bool`
        If `True` the :class:`~pyhdx.pyhdx.KineticSeries` instance is made uniform
    **metadata
        Dictionary of optional metadata. By default, holds the `temperature` and `pH` parameters.


    Attributes
    ----------
    state : :obj:`str`
        State of the kinetic series
    timepoints : :class:`~numpy.ndarray`
        Array with exposure times (sorted)
    peptides: :obj:`list`
        List of :class:`~pyhdx.models.PeptideMeasurements`, one list element per timepoint.
    cov: :class:`~pyhdx.models.Coverage`
        Coverage object describing peptide layout. When this `uniform` is `False`, this attribute is `None`

    """
    def __init__(self, data, **metadata):
        self.metadata = metadata
        assert len(np.unique(data['state'])) == 1
        self.state = data['state'][0]
        self.timepoints = np.sort(np.unique(data['exposure']))

        data_list = [(data[data['exposure'] == exposure]) for exposure in self.timepoints]
        sets = [{tuple(elem) for elem in fields_view(d, ['_start', '_end'])} for d in data_list]
        intersection = set.intersection(*sets)
        intersection_array = np.array([tup for tup in intersection], dtype=[('_start', int), ('_end', int)])

        # Select entries in data array which are in the interesection between all timepoints
        selected = [elem[np.isin(fields_view(elem, ['_start', '_end']), intersection_array)] for elem in data_list]
        self.peptides = [PeptideMeasurements(elem) for elem in selected]

        # Create coverage object from the first time point (as all are now equal)
        self.cov = Coverage(selected[0], c_term=metadata.get('c_term'))

    @property
    def temperature(self):
        try:
            return self.metadata['temperature']
        except KeyError:
            return None

    @temperature.setter
    def temperature(self, value):
        self.metadata['temperature'] = value

    @property
    def pH(self):
        try:
            return self.metadata['pH']
        except KeyError:
            return None

    @pH.setter
    def pH(self, value):
        self.metadata['pH'] = value

    @property
    def c_term(self):
        return self.cov.protein.c_term

    @c_term.setter
    def c_term(self, value):
        raise NotImplementedError('Cannot change c_term after object initialization')
        #todo allow this by making an new protein / coverage object
        self.metadata['c_term'] = value

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
    def scores_stack(self):
        """uptake scores to fit in a 2d stack"""
        scores_2d = np.stack([v.scores_average for v in self])
        return scores_2d

    @property
    def scores_peptides(self):
        scores_peptides = np.stack([v.scores for v in self])
        return scores_peptides

    @property
    def uptake_corrected(self):
        """matrix shape  N_t, N_p"""
        uptake_corrected = np.stack([v.uptake_corrected for v in self])
        return uptake_corrected


class PeptideMeasurements(Coverage):
    """
    Class with subset of peptides corresponding to only one state and exposure

    Parameters
    ----------
    data : :class`~numpy.ndarray`
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

    def __init__(self, data):
        assert len(np.unique(data['exposure'])) == 1, 'Exposure entries are not unique'
        assert len(np.unique(data['state'])) == 1, 'State entries are not unique'

        super(PeptideMeasurements, self).__init__(data)

        self.state = self.data['state'][0]
        self.exposure = self.data['exposure'][0]

    @property
    def scores(self):
        try:
            return self.data['scores']
        except ValueError:
            return self.data['uptake']

    @property
    def uptake_corrected(self):
        return self.data['uptake_corrected']

    @property
    def name(self):
        return self.state + '_' + str(self.exposure)

    @property
    def scores_average(self):
        return self.Z_norm.T.dot(self.scores)

    def calc_scores(self, residue_scores):
        """
        Calculates uptake scores per peptide given an array of individual residue scores

        Parameters
        ----------
        residue_scores : :class:`~numpy.ndarray`
            Array of scores per residue of length `prot_len`

        Returns
        -------

        scores : :class`~numpy.ndarray`
            Array of scores per peptide
        """

        scores = self.Z.dot(residue_scores)
        return scores


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


def series_intersection(series_list, fields=None):
    """
    Finds the intersection between peptides in :class:`~pydhx.models.KineticSeries` and returns new objects such that
    all peptides (coverage, exposure) between the series is identical.

    Optionally intersections by custom fields can be made.

    Parameters
    ----------
    series_list: :obj:`list`
        Input list of :class:`~pyhdx.models.KineticSeries`
    fields: :obj:`list`
        By which fields to take the intersections. Default is ['_start', '_end', 'exposure']

    Returns
    -------
    series_out: :obj:`list`
        Output list of :class:`~pyhdx.models.KineticSeries`
    """

    fields = fields or ['_start', '_end', 'exposure']

    full_arrays = [series.full_data for series in series_list]
    intersection = reduce(np.intersect1d, [fields_view(d, fields) for d in full_arrays])
    selected = [elem[np.isin(fields_view(elem, fields), intersection)] for elem in full_arrays]

    series_out = [KineticsSeries(data, **series.metadata) for data, series in zip(selected, series_list)]
    return series_out
