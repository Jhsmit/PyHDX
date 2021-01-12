from pyhdx.fileIO import read_dynamx, txt_to_protein
from pyhdx import PeptideMasterTable, KineticsFitting
import numpy as np
import pickle
from pathlib import Path
import torch

torch.manual_seed(43)
np.random.seed(43)
epochs = 1000

directory = Path(__file__).parent
fpath = directory / 'test_data' / 'ecSecB_apo.csv'

guess = True
state = 'SecB WT apo'
control = ('Full deuteration control', 0.167)

data = read_dynamx(fpath)

pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pf.set_control(control)
states = pf.groupby_state()
series = states[state]

temperature, pH = 273.15 + 30, 8.
kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

if guess:
    wt_avg_result = kf.weighted_avg_fit()
    output = wt_avg_result.output
    output.to_file(directory / 'test_data' / 'ecSecB_guess.txt')
else:
    output = txt_to_protein(directory / 'test_data' / 'ecSecB_guess.txt')


fr_torch = kf.global_fit_torch(output, epochs=epochs)
fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.txt')
with open(directory / 'test_data' / 'ecSecB_torch_fit.pick', 'wb') as f:
    pickle.dump(fr_torch, f)

