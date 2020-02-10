# -*- coding: utf-8 -*-

"""Main module."""

import numpy as np
import itertools
import scipy
from .math import solve_nnls

CSV_DTYPE = [
    ('Protein', 'U', 'c'),
    ('start', 'i', 'i'),
    ('end', 'i', 'i'),
    ('sequence', 'U', 'c'),
    ('modification', 'U', 'c'),
    ('fragment', 'U', 'c'),
    ('max_uptake', 'i', 'i'),
    ('MHP', 'f', 'f'),
    ('state', 'U', 'c'),
    ('exposure', 'f', 'f'),
    ('center', 'f', 'f'),
    ('center_sd', 'f', 'f'),
    ('uptake', 'f', 'f'),
    ('uptake_sd', 'f', 'f'),
    ('RT', 'f', 'f'),
    ('RT_sd', 'f', 'f')
]

HEADER = 'Protein,Start,End,Sequence,Modification,Fragment,MaxUptake,MHP,State,Exposure,Center,Center SD,Uptake,Uptake SD,RT,RT SD'


class PeptideCSVFile(object):
    """
    Input object of DynamX HDX .csv file

    Parameters
    ----------
    file_path : :obj:`str`
        File path of the .csv file
    drop_first : :obj:`int`
        Number of N-terminal amino acids to ignore. Default is 1.
    """
    def __init__(self, file_path, drop_first=1):

        names = [t[0] for t in CSV_DTYPE]
        self.data = np.genfromtxt(file_path, skip_header=1, delimiter=',', dtype=None, names=names, encoding='UTF-8')
        if drop_first:
            self.data['start'] += drop_first
            size = np.max([len(s) for s in self.data['sequence']])
            new_seq = np.array([s[drop_first:] for s in self.data['sequence']], dtype='<U' + str(size - drop_first))
            self.data['sequence'] = new_seq

    def return_by_name(self, control_state, control_exposure):
        #todo return dictionary of kinetic series instead

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
                out[name] = PeptideMeasurements(data_final, score)

        return out

    def get_data(self, state, exposure):
        """
        Get all data matching the supplied state and exposure.

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


class PeptideMeasurements(object):
    """
    Class with subset of peptides corresponding to only one state and exposure

    Parameters
    ----------
    data : :class`~numpy.ndarray`
        Numpy structured array with in put data
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

    def __init__(self, data, scores=None):
        assert len(np.unique(data['exposure'])) == 1, 'Exposure entries are not unique'
        assert len(np.unique(self.data['state'])) == 1, 'State entries are not unique'
        if scores is not None:
            assert len(scores) == len(data), 'Length of scores must match the length of the data (number of peptides)'
        self.data = data
        self.start = np.min(self.data['start'])
        self.stop = np.max(self.data['end'])  #todo refactor to end
        self.prot_len = self.stop - self.start + 1  # Total number of amino acids described by these measurments

        self.exposure = self.data['exposure'][0]
        self.state = self.data['state'][0]
        self.scores = scores #if scores is not None else None  # TODO currently needs scores for len

        self.big_X = np.zeros((len(data), self.prot_len), dtype=float)
        for row, entry in enumerate(self.data):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start
            pep_len = i1 - i0 + 1
            assert len(entry['sequence']) == pep_len
            self.big_X[row][i0:i1 + 1] = 1/pep_len

        ## sectioning
        combinations, num_letters = self.get_combinations(len(data))
        uid_mx = np.zeros((len(data), self.prot_len), dtype='U' + str(num_letters + 4))
        for c, (row, entry) in zip(combinations, enumerate(data)):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start
            # pep_len = i1 - i0 + 1
            uid_mx[row][i0:i1 + 1] = c

        big_sum = np.array([''.join(row) for row in uid_mx.T])
        b = big_sum == '' # Booleans where there is a gap
        regions = contiguous_regions(b)
        if len(regions) > 0:
            combinations, n = self.get_combinations(len(regions), prefix='gap_')
            for r, c in zip(regions, combinations):
                big_sum[r[0]:r[1]] = c

        vals, idx, counts = np.unique(big_sum, return_index=True, return_counts=True)
        vals = vals[np.argsort(idx)]  #
        self.counts = counts[np.argsort(idx)]
        self.cs = np.cumsum(self.counts)

        #states is one if coverage, 0 if not
        self.states = np.ones_like(self.counts)
        gaps = vals.astype('<U4') == 'gap_'
        self.states[gaps] = 0

        # Matrix of N columns with N equal to sections of identical residues.
        # Values are the number of residues in the blocks
        self.X = np.zeros((len(data), len(self.counts)), dtype=float)
        for row, entry in enumerate(data):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start + 1
            p = np.searchsorted(self.cs, [i0, i1], side='right')

            self.X[row][p[0]:p[1]] = self.counts[p[0]:p[1]]

    def __len__(self):
        return len(self.data)

    @staticmethod
    def get_combinations(num, prefix=''):
        """returns unique combinations of letters"""
        abc = 'abcdefghijklmnopqrstuvwxyz'
        num_letters = int(np.ceil(np.log(num) / np.log(len(abc))))
        combinations = list([prefix + ''.join(letters) for letters in itertools.combinations(abc, num_letters)])

        return combinations[:num], num_letters

    @property
    def name(self):
        return self.state + '_' + str(self.exposure)

    @property
    def X_norm(self):
        return self.X / np.sum(self.X, axis=1)[:, np.newaxis]

    @property
    def big_X_norm(self):
        return self.big_X / np.sum(self.big_X, axis=0)[np.newaxis, :]

    @property
    def big_X_sq_norm(self):
        return self.big_X**2 / np.sum(self.big_X**2, axis=0)[np.newaxis, :]

    @property
    def scores_average(self):
        return self.big_X_norm.T.dot(self.scores)

    @property
    def scores_average_sq(self):
        return self.big_X_sq_norm.T.dot(self.scores)

    @property
    def scores_lstsq(self):
        x, res, rank, s = np.linalg.lstsq(self.X_norm, self.scores)
        return np.repeat(x, self.counts)

    def scores_nnls_tikonov(self, reg):
        x = solve_nnls(self.X_norm.T, self.scores, reg=reg)
        return np.repeat(x, self.counts)

    def scores_nnls(self):
        x = scipy.optimize.nnls(self.X_norm, self.scores,)[0]
        return np.repeat(x, self.counts)

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

        scores = self.big_X.dot(residue_scores)
        return scores

    @property
    def sequence(self):
        """:obj:`str: String of the full protein sequence. Gaps of no coverage are filled with Prolines."""
        seq = np.full(self.stop, 'P', dtype='U')
        for d in self.data:
            i = d['start'] - 1
            j = d['end']
            seq[i:j] = [s for s in d['sequence']]
        return ''.join(seq)


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
