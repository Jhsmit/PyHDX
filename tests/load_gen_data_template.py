import os
import numpy as np
from pyhdx.support import np_from_txt, fmt_export
from pyhdx import PeptideMasterTable, KineticsFitting
import pickle

fit_dir = 'test_data'
directory = os.path.dirname(__file__)
np.random.seed(43)


fpath = os.path.join(directory, 'test_data', 'simulated_data.csv')
data = np_from_txt(fpath, delimiter=',')
data['end'] += 1  # because this simulated data is in old format of inclusive, inclusive
sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
start, end = 5, 45  # total span of protein (inc, inc)
nc_start, nc_end = 31, 34  # span of no coverage area (inc, inc)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
states = pmt.groupby_state()
series = states['state1']

print(series.scores_peptides.T.shape)

print(series.cov.X.shape)

init_arr = np_from_txt(os.path.join(fit_dir, 'fit_simulated_wt_avg.txt'))
