import numpy as np
from numpy.lib.recfunctions import append_fields
import itertools
import scipy
from functools import reduce
from operator import add
from pyhdx.math import solve_nnls
from pyhdx.support import reduce_inter, make_view
from pyhdx.expfact.kint import calculate_kint_per_residue


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
    ignore_prolines: :obj:`bool`
        Boolean to toggle ignoring of proline residues. When True these residues are treated as if they're not present
        in the protein.
    sort: :obj:`bool`
        Set to ``True`` to sort the input. Sort order is 'start', 'end', 'sequence', 'exposure', 'state'.
    remove_nan: :obj`bool`
        Set to ``True`` to remove NaN entries in uptake

    """

    def __init__(self, data, drop_first=1, ignore_prolines=True, sort=True, remove_nan=True):
        assert np.all(data['start'] < data['end']), 'All `start` entries must be smaller than their `end` entries'
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


        #Covert sequence to upper case if not so already
        self.data['sequence'] = [s.upper() for s in self.data['sequence']]
        #Mark ignored prolines with lower case letters
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

        ex_residues = [len(s) - s.count('x') - s.count('p') for s in self.data['sequence']]
        if 'ex_residues' not in self.data.dtype.names:
            self.data = append_fields(self.data, ['ex_residues'], [ex_residues], usemask=False)

    def __len__(self):
        return len(self.data)

    def groupby_state(self, make_uniform=True):
        """
        Groups measurements in the dataset by state and returns them in a dictionary as a
        :class:`~pyhdx.pyhdx.KineticSeries`.

        Parameters
        ----------
        make_uniform : :obj:`bool`
            If `True` the returned :class:`~pyhdx.pyhdx.KineticSeries` is made uniform

        Returns
        -------
        out : :obj:`dict`
            Dictionary where keys are state names and values are :class:`~pyhdx.pyhdx.KineticSeries`
        """

        states = np.unique(self.data['state'])
        return {state: KineticsSeries(self.data[self.data['state'] == state], make_uniform=make_uniform) for state in states}

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

        #raise NotImplementedError()

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


class TFCoverage(object):
    """
    Object describing layout and coverage of peptides and generating the corresponding matrices. Peptides should all
    belong to the same state and have the same exposure time.

    Parameters
    ----------
    data : ~class:`~numpy.ndarray`
        Numpy structured array with input peptides

    Attributes
    ----------
        start : :obj:`int`
            Index of residue first appearing in the peptides (first residue is 1)
        end : :obj:`int`
            Index of last residue appearing in the peptides (inclusive)
        r_number : :class:`~numpy.ndarray`
            Array of residue numbers which are covered by the peptides, excluding prolines if they are set to be ignored.
        prot_len : :obj:`int`
            Total number of residues the peptides covering, excluding prolines if they are set to be ignored.
        X : :class:`~numpy.ndarray`
            N x M matrix where N is the number of peptides and M equal to `prot_len`.
            Values are 1/(ex_residues) where there is coverage, so that rows sum to 1
    """

    def __init__(self, data):
        assert len(np.unique(data['exposure'])) == 1, 'Exposure entries are not unique'
        assert len(np.unique(data['state'])) == 1, 'State entries are not unique'


        self.data = data
        self.start = np.min(self.data['start'])
        self._start = np.min(self.data['_start'])
        self.end = np.max(self.data['end'])
        self._end = np.max(self.data['_end'])
        self.r_number = np.arange(self.start, self.end)

        self.prot_len = len(self.r_number)

        self.X = np.zeros((len(self.data), len(self.r_number)), dtype=np.int)  # cast to float for tf fit
        for row, entry in enumerate(self.data):
            i0, i1 = np.searchsorted(self.r_number, (entry['start'], entry['end']))
            self.X[row][i0:i1] = 1

    @property
    def X_norm(self):
        """:class:`~np.ndarray`: `X` coefficient matrix normalized column wise."""
        return self.X / np.sum(self.X, axis=0)[np.newaxis, :]

    @property
    def has_coverage(self):
        """:class:`~np.ndarray`: Boolean array indicating if the residues along r_number have coverage"""
        return np.sum(self.X, axis=0) > 0

    def __len__(self):
        return len(self.data)

    @property  #todo refactor to sequence
    def cov_sequence(self):
        """amino acids one letter codes corresponding to r_number array"""
        idx = np.searchsorted(self.sequence_r_number, self.r_number)
        seq_arr = np.array([s for s in self.sequence])
        return seq_arr[idx]

    @property  #todo refactor to full sequence. Perhaps also return as array
    def sequence(self):
        """:obj:`str`: String of the full protein sequence. One letter coding where X marks regions of no coverage"""

        r_number = self.sequence_r_number
        seq = np.full_like(r_number, fill_value='X', dtype='U')

        for d in self.data:
            i, j = np.searchsorted(r_number, [d['_start'], d['_end']])
            seq[i:j] = [s for s in d['_sequence']]
        return ''.join(seq)

    @property
    def sequence_r_number(self):  #todo refactor to: full_r_number
        """~class:`numpy.ndarray`: Array of r numbers corresponding to residues in sequence"""
        start = min(self._start, 1)  # start at least at 1 unless the protein extends into negative numbers
        r_number = np.arange(start, self._end)

        return r_number

    def calc_kint(self, temperature, pH, c_term):
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
        c_term : :obj:`int`
            index of the last residue in the sequence (first residue is 1)

        Returns
        -------

        k_int : ~class:`~numpy.ndarray`
            Array of intrisic exchange rates

        """

        c_term = len(self.sequence) + 1 if c_term is None else c_term
        k_int_list = [-1.]
        for i, (previous, current) in enumerate(zip(self.sequence[:-1], self.sequence[1:])):
            if previous == 'X' or current == 'X':
                k_int_list.append(0.)
            elif current == 'P':
                k_int_list.append(0.)
            else:
                k_int = calculate_kint_per_residue(previous, current, i + 2, c_term, temperature, pH)
                k_int_list.append(k_int)
        return np.array(k_int_list)

    def get_kint_array(self, temperature, pH, c_term):
        k_int = self.calc_kint(temperature, pH, c_term=c_term)
        k_r_number = self.sequence_r_number

        output = np.empty(len(k_int), dtype=[('r_number', np.int), ('k_int', np.float)])
        output['r_number'] = k_r_number
        output['k_int'] = k_int

        return output

    def split(self, gap_size=-1):
        """
        Splits the dataset into independent parts which have no overlapping peptides between them. To determine overlap,
        the modified 'start' and 'end' fields are used which take into account N-terminal non-exchanging residues and
        prolines.

        Returns
        -------

        output: :obj:`dict`
            Dictionary where keys are {start}_{end} (inclusive, exclusive) of the corresponding sections, values are
            of the ``type`` of the current instance.
        gap_size: :obj:`int`
            Gaps of this size between adjacent peptides is not considered to overlap. A value of -1 means that peptides
            with exactly zero overlap are separated. With gap_size=0 peptides with exactly zero overlap are not separated,
            and larger values tolerate larger gap sizes.
        """

        klass = self.__class__

        # intervals = [(s, e + 1) for s, e in zip(self.data['start'], self.data['end'])]
        # sections = reduce_inter(intervals, gap_size=gap_size)

        sections = self.get_sections(gap_size)

        output = {}
        for s, e in sections:

            b = np.logical_and(self.data['start'] >= s, self.data['end'] <= e)
            output[f'{s}_{e}'] = klass(self.data[b])

        return output

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
        start : :obj:`int`
            Index of residue first appearing in the peptides (first residue is 1)
        end : :obj:`int`
            Index of last residue appearing in the peptides (inclusive)
        r_number : :class:`~numpy.ndarray`
            Array of residue numbers which are covered by the peptides, excluding prolines if they are set to be ignored.
        prot_len : :obj:`int`
            Total number of residues the peptides covering, excluding prolines if they are set to be ignored.
        X : :class:`~numpy.ndarray`
            N x M matrix where N is the number of peptides and M equal to `prot_len`.
            Values are 1/(ex_residues) where there is coverage, so that rows sum to 1
    """

    def __init__(self, data):
        assert len(np.unique(data['exposure'])) == 1, 'Exposure entries are not unique'
        assert len(np.unique(data['state'])) == 1, 'State entries are not unique'
        # character to use to count for proline occurences

        # todo insert and update coverage logic

        self.data = data
        self.start = np.min(self.data['start'])
        self._start = np.min(self.data['_start'])
        self.end = np.max(self.data['end'])
        self._end = np.max(self.data['_end'])
        self.r_number = np.arange(self.start, self.end)

        # Find all indices of prolines in the middle of sequences, remove from r_number array and from sequence
        p = [entry['_start'] + i - self.start for entry in self.data for i, s in enumerate(entry['sequence']) if s == 'p']
        p_index = np.unique(p)
        self.r_number = np.delete(self.r_number, p_index)  # remove r number indices
        self.prot_len = len(self.r_number)

        self.X = np.zeros((len(self.data), len(self.r_number)), dtype=float)
        for row, entry in enumerate(self.data):
            i0, i1 = np.searchsorted(self.r_number, (entry['start'], entry['end']))
            self.X[row][i0:i1] = 1 / entry['ex_residues']

    def __len__(self):
        return len(self.data)

    @property
    def block_length(self):
        """:class:`~numpy.ndarary`: Lengths of unique blocks of residues in the peptides map,
            along the `r_number` axis"""

        # indices are start and stop values of blocks
        print(self.data['end'])
        indices = np.sort(np.concatenate([self.data['start'], self.data['end']]))
        print(indices)
        #indices of insertion into r_number vector gives us blocks with taking prolines into account.
        diffs = np.diff(np.searchsorted(self.r_number, indices))

        block_length = diffs[diffs != 0]
        return block_length

    @property
    def block_coverage(self):
        """:class:`~np.ndarray`: Boolean array with ``True`` values where blocks have coverage."""
        block_coverage = np.sum(self.X_red, axis=0) > 0
        return block_coverage

    @property
    def X_red(self):
        """:class:`~np.ndarray`: Reduced NxM coefficient matrix, with N the number of peptides and M the number of blocks.
            Elements are equal to the length of the block.
        """
        cs = np.cumsum(self.block_length)
        X_red = np.zeros((len(self.data), len(self.block_length)), dtype=float)
        for row, entry in enumerate(self.data):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start
            p = np.searchsorted(cs, [i0, i1], side='right')

            X_red[row][p[0]:p[1]] = self.block_length[p[0]:p[1]]

        return X_red

    @property
    def X_red_norm(self):
        """:class:`~np.ndarray`: `X_red` blocks coefficient matrix normalized column wise."""
        return self.X_red / np.sum(self.X_red, axis=0)[np.newaxis, :]

    @property
    def X_norm(self):
        """:class:`~np.ndarray`: `X` coefficient matrix normalized column wise."""
        return self.X / np.sum(self.X, axis=0)[np.newaxis, :]

    @property
    def has_coverage(self):
        """:class:`~np.ndarray`: Boolean array indicating if the residues along r_number have coverage"""
        return np.sum(self.X, axis=0) > 0



    @property
    def sequence(self):
        """:obj:`str`: String of the full protein sequence. One letter coding where X marks regions of no coverage"""

        r_number = self.sequence_r_number
        seq = np.full_like(r_number, fill_value='X', dtype='U')

        for d in self.data:
            i, j = np.searchsorted(r_number, [d['_start'], d['_end']])
            seq[i:j] = [s for s in d['_sequence']]
        return ''.join(seq)

    @property
    def sequence_r_number(self):
        """~class:`numpy.ndarray`: Array of r numbers corresponding to residues in sequence"""
        start = min(self._start, 1)  # start at least at 1 unless the protein extends into negative numbers
        r_number = np.arange(start, self._end)

        return r_number

    def calc_kint(self, temperature, pH, c_term):
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
        c_term : :obj:`int`
            index of the last residue in the sequence (first residue is 1)

        Returns
        -------

        k_int : ~class:`~numpy.ndarray`
            Array of intrisic exchange rates

        """

        c_term = len(self.sequence) + 1 if c_term is None else c_term
        k_int_list = [-1.]
        for i, (previous, current) in enumerate(zip(self.sequence[:-1], self.sequence[1:])):
            if previous == 'X' or current == 'X':
                k_int_list.append(0.)
            elif current == 'P':
                k_int_list.append(-1)
            else:
                k_int = calculate_kint_per_residue(previous, current, i + 2, c_term, temperature, pH)
                k_int_list.append(k_int)
        return np.array(k_int_list)

    def split(self, gap_size=-1):
        """
        Splits the dataset into independent parts which have no overlapping peptides between them. To determine overlap,
        the modified 'start' and 'end' fields are used which take into account N-terminal non-exchanging residues and
        prolines.

        Returns
        -------

        output: :obj:`dict`
            Dictionary where keys are {start}_{end} (inclusive, exclusive) of the corresponding sections, values are
            of the ``type`` of the current instance.
        gap_size: :obj:`int`
            Gaps of this size between adjacent peptides is not considered to overlap. A value of -1 means that peptides
            with exactly zero overlap are separated. With gap_size=0 peptides with exactly zero overlap are not separated,
            and larger values tolerate larger gap sizes.
        """

        klass = self.__class__

        # intervals = [(s, e + 1) for s, e in zip(self.data['start'], self.data['end'])]
        # sections = reduce_inter(intervals, gap_size=gap_size)

        sections = self.get_sections(gap_size)

        output = {}
        for s, e in sections:

            b = np.logical_and(self.data['start'] >= s, self.data['end'] <= e)
            output[f'{s}_{e}'] = klass(self.data[b])

        return output

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

    Attributes
    ----------
    state : :obj:`str`
        State of the kinetic series
    timepoints : :class:`~numpy.ndarray`
        Array with exposure times (sorted)

    """
    def __init__(self, data, make_uniform=True, **metadata):
        if isinstance(data, np.ndarray):
            assert len(np.unique(data['state'])) == 1
            self.state = data['state'][0]
            self.timepoints = np.sort(np.unique(data['exposure']))

            self.peptides = [PeptideMeasurements(data[data['exposure'] == exposure]) for exposure in self.timepoints]

            if make_uniform:
                self.make_uniform()

            if self.uniform:
                self.cov = Coverage(data[data['exposure'] == self.timepoints[0]])
            else:
                self.cov = None

        elif isinstance(data, list):
            for elem in data:
                assert isinstance(elem, PeptideMeasurements), 'Invalid data format'
                assert elem.state == data[0].state, 'All peptide states must be equal'
            data_sort = sorted(data, key=lambda x: x.exposure)
            self.state = data[0].state
            self.timepoints = np.array([elem.exposure for elem in data_sort])
            self.peptides = data

            if self.uniform:
                # Use first peptideset to create coverage object
                self.cov = Coverage(data[0].data)
                self.tf_cov = TFCoverage(data[0].data)
            else:
                self.cov = None

        else:
            raise TypeError('Invalid data type')

    def make_uniform(self, in_place=True):
        """
        Removes entries from time points, ensuring that all time points have equal coverage


        Returns
        -------

        """
        #todo perhaps move to a function

        sets = [{(s, e, seq) for s, e, seq in zip(pm.data['start'], pm.data['end'], pm.data['sequence'])} for pm in self]
        intersection = set.intersection(*sets)
        dtype = [('start', int), ('end', int), ('sequence', self.full_data['sequence'].dtype)]
        inter_arr = np.array([tup for tup in intersection], dtype=dtype)

        if in_place:
            self.peptides = [pm[np.isin(pm.data[['start', 'end', 'sequence']], inter_arr)] for pm in self]
            self.cov = Coverage(self[0].data)
            self.tf_cov = TFCoverage(self[0].data)


        else:
            raise NotImplementedError('Only making peptidesets uniform in place is implemented')

    @property
    def uniform(self):
        """Returns ``True`` if for all time point coverages are equal"""
        is_uniform = np.all([self[0] == elem for elem in self])
        return is_uniform

    @property
    def full_data(self):
        """returns the full dataset of all timepoints"""
        full_data = np.concatenate([pm.data for pm in self])
        return full_data

    def split(self, gap_size=-1):
        """
        Splits the dataset into independent parts which have no overlapping peptides between them

        Returns
        -------
        output : :obj:`dict`
            Output dictionary with individual kinetic series. Keys are '{start}_{stop}', (including, excluding)
             values are :class:`~pyhdx.pyhdx.KineticSeries` objects.
        gap_size: :obj:`int`
            Gaps of this size between adjacent peptides is not considered to overlap. A value of -1 means that peptides
            with exactly zero overlap are separated. With gap_size=0 peptides with exactly zero overlap are not separated,
            and larger values tolerate larger gap sizes.

        """
        if self.uniform:
            # end is inclusive therefore +1 is needed
            intervals = [(s, e + 1) for s, e in zip(self[0].data['start'], self[0].data['end'])]
        else:
            raise AssertionError("not uniform data not yet supported")
            intervals = reduce(add, [[(s, e) for s, e in zip(pm.data['start'], pm.data['end'])] for pm in self])


        split_list = [pm.split(gap_size=gap_size) for pm in self]
        #accumulate all keys in the split list and sort them by start then end
        keys = sorted(np.unique([list(dic.keys()) for dic in split_list]), key=lambda x: tuple(int(c) for c in x.split('_')))
        #keys = ''
        #sections = reduce_inter(intervals)
        output = {}
        for key in keys:
            output[key] = KineticsSeries(list([dic[key] for dic in split_list]))
            # s, e = section
            # b = np.logical_and(full_ds['start'] >= s, full_ds['end'] <= e)
            # output['{}_{}'.format(s, e)] = KineticsSeries(full_ds[b])

        return output

    def __len__(self):
        return len(self.timepoints)

    def __iter__(self):
        return self.peptides.__iter__()

    def __getitem__(self, item):
        return self.peptides.__getitem__(item)

    def set_control(self, control_100, control_zero=None, remove_nan=True):
        """
        Apply a control dataset to the underlying PeptideMeasurements of this object. A `scores` attribute is added to
        the PeptideMeasurement by normalizing its uptake value with respect to the control uptake value to 100%. Entires
        which are in the measurement and not in the control or vice versa are deleted.
        Optionally, ``control_zero`` can be specified which is a datasets whose uptake value will be set to zero.

        Parameters
        ----------
        control_100 : :class:`~numpy.ndarray`
            Numpy structured array with control peptides to use for normalization to 100%
        control_zero : :class:`~numpy.ndarray`
            Numpy structured array with control peptides to use for normalization to 0%
        remove_nan : :obj:`Bool`
            If `True`, `NaN` entries are removed from the controls
        Returns
        -------

        """

        raise DeprecationWarning('will be removed')

        for pm in self:
            pm.set_control(control_100, control_zero, remove_nan=remove_nan)

    @property
    def scores_stack(self):
        # todo move this to series
        """uptake scores to fit in a 2d stack"""
        scores_2d = np.stack([v.scores_average for v in self])
        return scores_2d

    @property
    def scores_norm(self):
        # Normalized to 100 array of scores
        print('where is this used?')
        scores_norm = 100 * (self.scores_stack / self.scores_stack[-1, :][np.newaxis, :])
        return scores_norm

    @property
    def scores_peptides(self):
        scores_peptides = np.stack([v.scores for v in self])
        return scores_peptides

    @property
    def uptake_corrected(self):
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
      #  self.scores = data['uptake']

    @property
    def scores(self):
        try:
            return self.data['scores']
        except ValueError:
            return self.data['uptake']

    @property
    def uptake_corrected(self):
        return self.data['uptake_corrected']

    def set_control(self, control_100, control_0=None, remove_nan=True):
        """
        Apply a control dataset to this object. A `scores` attribute is added to the object by normalizing its uptake
        value with respect to the control uptake value to 100%. Entries which are in the measurement and not in the
        control or vice versa are deleted.
        Optionally, ``control_zero`` can be specified which is a datasets whose uptake value will be set to zero.

        Parameters
        ----------
        control_100 : :class:`~numpy.ndarray`
            Numpy structured array with control peptides to use for normalization to 100%
        control_0 : :class:`~numpy.ndarray`
            Numpy structured array with control peptides to use for normalization to 0%
        remove_nan : :obj:`Bool`
            If `True`, `NaN` entries are removed from the controls

        Returns
        -------

        """
        # peptides in measurements that are also in the control
        #todo check for NaNs with lilys file
        raise DeprecationWarning("will be deprecated")

        if control_0 is None:
            control_0 = np.copy(control_100)
            control_0['uptake'] = 0

        # Remove NaN entries from controls
        if remove_nan:
            control_100 = control_100[~np.isnan(control_100['uptake'])]
            control_0 = control_0[~np.isnan(control_0['uptake'])]

        b_100 = np.isin(self.data['sequence'], control_100['sequence'])
        b_0 = np.isin(self.data['sequence'], control_0['sequence'])
        data_selected = self.data[np.logical_and(b_100, b_0)]

        # Control peptides corresponding to those peptides in measurement
        c_100_selected = control_100[np.isin(control_100['sequence'], data_selected['sequence'])]
        c_0_selected = control_0[np.isin(control_0['sequence'], data_selected['sequence'])]

        # Sort both datasets by starting index and then by sequence to make sure they are both equal
        data_final = np.sort(data_selected, order=['start', 'sequence'])
        control_100_final = np.sort(c_100_selected, order=['start', 'sequence'])
        control_0_final = np.sort(c_0_selected, order=['start', 'sequence'])

        #todo move assert to testing
        assert np.all(data_final['sequence'] == control_100_final['sequence'])
        assert np.all(data_final['start'] == control_100_final['start'])
        assert np.all(data_final['end'] == control_100_final['end'])

        scores = 100 * ( (data_final['uptake'] - control_0_final['uptake']) /
                (control_100_final['uptake'] - control_0_final['uptake']) )

        #update this when changing to Coverage objects

        super(PeptideMeasurements, self).__init__(data_final)
        self.scores = scores

    def __getitem__(self, item):
        if isinstance(item, int):
            return None
        else:
            data = self.data[item]
            scores = self.scores[item]

            pm = PeptideMeasurements(data)
#            pm.scores = scores
            return pm

    @property
    def name(self):
        return self.state + '_' + str(self.exposure)

    @property
    def scores_average(self):
        return self.X_norm.T.dot(self.scores)

    @property
    def scores_lstsq(self):
        """DEPRECATED"""
        x, res, rank, s = np.linalg.lstsq(self.X_norm, self.scores)
        return np.repeat(x, self.block_length)

    def scores_nnls_tikonov(self, reg):
        """DEPRECATED"""
        x = solve_nnls(self.X_norm.T, self.scores, reg=reg)
        return np.repeat(x, self.block_length)

    def scores_nnls(self):
        """DEPRECATED"""
        x = scipy.optimize.nnls(self.X_norm, self.scores,)[0]
        return np.repeat(x, self.block_length)

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

        scores = self.X.dot(residue_scores)
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
