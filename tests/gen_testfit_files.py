import os
import numpy as np
from pyhdx import PeptideMasterTable, KineticsFitting, read_dynamx
from pathlib import Path
import pickle
import torch

directory = Path(__file__).parent

torch.manual_seed(43)
np.random.seed(43)
epochs = 1000

fpath = directory / 'test_data' / 'simulated_data_uptake.csv'
data = read_dynamx(fpath)
sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
start, end = 5, 45  # total span of protein (inc, inc)
nc_start, nc_end = 31, 34  # span of no coverage area (inc, inc)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_backexchange(0.)
states = pmt.groupby_state()
series = states['state1']

temperature, pH = 300, 8
series.cov.protein.set_k_int(temperature=temperature, pH=pH)
series.cov.protein.to_file(directory / 'test_data' / 'simulated_data_info.txt')

kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

fr1 = kf.weighted_avg_fit()
out1 = fr1.output

out1.to_file(directory / 'test_data' / 'fit_simulated_wt_avg.txt')
with open(directory / 'test_data' / 'fit_simulated_wt_avg.pick', 'wb') as f:
    pickle.dump(fr1, f)

fr_torch = kf.global_fit_torch(out1, epochs=epochs)
out_deltaG = fr_torch.output

out_deltaG.to_file(directory / 'test_data' / 'fit_simulated_torch.txt')
with open(os.path.join(directory / 'test_data' / 'fit_simulated_torch.pick'), 'wb') as f:
    pickle.dump(fr_torch, f)
