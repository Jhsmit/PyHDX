import os
import numpy as np
from pyhdx import PeptideMasterTable, KineticsFitting, read_dynamx
from pyhdx.fileIO import txt_to_np
import pickle

fit_dir = 'test_data'
directory = os.path.dirname(__file__)
np.random.seed(43)

fpath = os.path.join(directory, 'test_data', 'ds1.csv')
data = read_dynamx(fpath)
sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
start, end = 5, 45  # total span of protein (inc, inc)
nc_start, nc_end = 31, 34  # span of no coverage area (inc, inc)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_backexchange(0.)
states = pmt.groupby_state()

#series = states['state1']
series = states['PpiANative']
series.make_uniform()

print(series.scores_peptides.T.shape)
print(series.uptake_corrected.shape)  ## N_t, N_p
print(series.cov.X.shape)

init_arr = txt_to_np(os.path.join(fit_dir, 'fit_simulated_wt_avg.txt'))
