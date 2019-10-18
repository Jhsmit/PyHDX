# -*- coding: utf-8 -*-

"""Main module."""

import numpy as np
import itertools
import scipy
from .math import solve_nnls

CSV_DTYPE = [
    ('Protein', 'U2'),
    ('start', 'i'),
    ('end', 'i'),
    ('sequence', 'U64'),
    ('modification', 'U12'),
    ('fragment', 'U'),
    ('max_uptake', 'i'),
    ('MHP', 'f'),
    ('state', 'U24'),
    ('exposure', 'f'),
    ('center', 'f'),
    ('center_sd', 'f'),
    ('uptake', 'f'),
    ('uptake_sd', 'f'),
    ('RT', 'f'),
    ('RT_sd', 'f')
]


class PeptideCSVFile(object):
    def __init__(self, file_path):
        self.data = np.genfromtxt(file_path, skip_header=1, delimiter=',', dtype=CSV_DTYPE)

    def return_by_name(self, control_state, control_exposure):
        bc1 = self.data['state'] == control_state
        bc2 = self.data['exposure'] == control_exposure

        control = self.data[np.logical_and(bc1, bc2)]
        st = np.unique(self.data['state'])
        exp = np.unique(self.data['exposure'])

        out = {}
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

            #TODO finalize name
            #bs = score < 100
            #out[name] = PeptideMeasurements(d[bs], score[bs])
            if len(score) > 0:
                out[name] = PeptideMeasurements(data_final, score)

        return out

    def return_measurements(self, i_control, i_max):
        # Deprecated
        others = [i for i in range(i_max) if i != i_control]

        control = self.data[i_control::i_max]['uptake']
        out = {}
        for i in others:
            d = self.data[i::i_max]
            score = 100*d['uptake'] / control
            name = d['state'][0] + '_' + str(d['exposure'][0])
            out[name] = PeptideMeasurements(d, score)

        return out


class PeptideMeasurements(object):
    def __init__(self, data, scores=None, min_len=0):
        """data: structured array with at least the fields 'start', 'end', 'sequence' 'exposure'

        """
        self.data = data
        self.start = np.min(self.data['start'])
        self.stop = np.max(self.data['end'])
        self.prot_len = self.stop - self.start + 1

        if len(np.unique(self.data['exposure'])) == 1:
            self.exposure = self.data['exposure'][0]
        else:
            self.exposure = None

        if len(np.unique(self.data['state'])) == 1:
            self.state = self.data['state'][0]
        else:
            self.state = None

        self.scores = scores if scores is not None else None  # TODO currently needs scores for len

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
        #big_sum[b] = combinations

        vals, idx, counts = np.unique(big_sum, return_index=True, return_counts=True)
        vals = vals[np.argsort(idx)]  #
        self.counts = counts[np.argsort(idx)]
        self.cs = np.cumsum(self.counts)

        #states is one if coverage, 0 if not
        self.states = np.ones_like(self.counts)
        gaps = vals.astype('<U4') == 'gap_'
        self.states[gaps] = 0

        # dont use this, doesnt work
        if min_len > 0:
            merged_counts = np.array([], dtype=int)
            rr_prev = 0
            for rl, rr in contiguous_regions(self.counts <= 2):
                merged_counts = np.append(merged_counts, self.counts[rr_prev:rl])
                rr_prev = rr
                s = self.counts[rl:rr].sum()
                merged_counts = np.append(merged_counts, s)

            final_counts = merged_counts.copy()
            for i, e in enumerate(merged_counts):
                if e <= 2:
                    final_counts[i] -= e
                    if i == 0:
                        final_counts[i + 1] += e
                    elif i == len(merged_counts) - 1:
                        final_counts[i - 1] += e
                    elif final_counts[i + 1] < final_counts[i - 1]:
                        final_counts[i + 1] += e
                    else:
                        final_counts[i - 1] += e
            final_counts = final_counts[final_counts > 0]
            assert final_counts.sum() == self.counts.sum()
            self.counts = final_counts
            cs = np.cumsum(self.counts)


        # Matrix of N columns with N equal to sections of identical residues.
        # Values are the number of residues in the blocks
        self.X = np.zeros((len(data), len(self.counts)), dtype=float)
        for row, entry in enumerate(data):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start + 1
            p = np.searchsorted(self.cs, [i0, i1], side='right')

            self.X[row][p[0]:p[1]] = self.counts[p[0]:p[1]]

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

    def get_scores(self, method):
        if method == 'avg':
            return self.scores_average

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
        print(np.any(np.isnan(self.X_norm)))
        print(np.any(np.isinf(self.X_norm)))
        x = scipy.optimize.nnls(self.X_norm, self.scores,)[0]
        return x

    def calc_scores(self, residue_scores):
        return self.big_X.dot(residue_scores)




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
