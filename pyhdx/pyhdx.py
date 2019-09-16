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
    ('state', 'U12'),
    ('exposure', 'f'),
    ('center', 'f'),
    ('center_sd', 'f'),
    ('uptake', 'f'),
    ('uptake_sd', 'f'),
    ('RT', 'f'),
    ('RT_sd', 'f')
]


class PeptideCSVFile(object):
    def __init__(self, file_path, i_control, i_max):
        self.data = np.genfromtxt(file_path, skip_header=1, delimiter=',', dtype=CSV_DTYPE)
        self.i_control = i_control
        self.i_max = i_max

    def return_measurements(self):
        others = [i for i in range(self.i_max) if i != self.i_control]

        control = self.data[self.i_control::self.i_max]['uptake']

        out = {}
        for i in others:
            d = self.data[i::self.i_max]
            score = 100*d['uptake'] / control
            name = d['state'][0] + '_' + str(d['exposure'][0])
            out[name] = PeptideMeasurements(d, score)

        return out


class PeptideMeasurements(object):
    def __init__(self, data, scores=None, min_len=0):
        """data: structured array with at least the fields 'start', 'end', 'sequence'

        """
        self.data = data
        self.start = np.min(self.data['start'])
        self.stop = np.max(self.data['end'])
        self.prot_len = self.stop - self.start + 1

        self.scores = scores if scores is not None else None  # TODO currently needs scores for len

        self.big_X = np.zeros((len(data), self.prot_len), dtype=float)
        for row, entry in enumerate(self.data):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start
            pep_len = i1 - i0 + 1
            assert len(entry['sequence']) == pep_len
            self.big_X[row][i0:i1 + 1] = 1/pep_len


        abc = 'abcdefghijklmnopqrstuvwxyz'
        num_letters = int(np.ceil(np.log(len(data)) / np.log(len(abc))))
        combinations = list([''.join(letters) for letters in itertools.combinations(abc, num_letters)])

        uid_mx = np.zeros((len(data), self.prot_len), dtype='U2')
        for c, (row, entry) in zip(combinations, enumerate(data)):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start
            pep_len = i1 - i0 + 1
            uid_mx[row][i0:i1 + 1] = c

        big_sum = np.array([''.join(row) for row in uid_mx.T])
        vals, idx, counts = np.unique(big_sum, return_index=True, return_counts=True)
        vals = vals[np.argsort(idx)]  #
        self.counts = counts[np.argsort(idx)]
        cs = np.cumsum(self.counts)

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
            print(final_counts)
            self.counts = final_counts
            cs = np.cumsum(self.counts)


        print(cs)
        # Matrix of N columns with N equal to sections of identical residues.
        # Values are the number of residues in the blocks
        self.X = np.zeros((len(data), len(self.counts)), dtype=float)
        for row, entry in enumerate(data):
            i0 = entry['start'] - self.start
            i1 = entry['end'] - self.start + 1
            p = np.searchsorted(cs, [i0, i1], side='right')
            print(row, p, i0, i1)

            self.X[row][p[0]:p[1]] = self.counts[p[0]:p[1]]

    @property
    def X_norm(self):
        return self.X / np.sum(self.X, axis=1)[:, np.newaxis]

    @property
    def big_X_norm(self):
        return self.big_X / np.sum(self.big_X, axis=0)[np.newaxis, :]

    def get_scores(self, method):
        if method == 'avg':
            return self.scores_average

    @property
    def scores_average(self):
        return self.big_X_norm.T.dot(self.scores)

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
