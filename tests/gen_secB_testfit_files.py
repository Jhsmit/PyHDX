from pyhdx.fileIO import read_dynamx, txt_to_protein
from pyhdx import PeptideMasterTable, KineticsFitting, BatchFitting
import numpy as np
import pickle
from pathlib import Path
import torch

torch.manual_seed(43)
np.random.seed(43)
epochs = 1000

directory = Path(__file__).parent
test_data_dir = directory / 'test_data'

guess = False
control = ('Full deuteration control', 0.167)

data = read_dynamx(test_data_dir / 'ecSecB_apo.csv', test_data_dir / 'ecSecB_dimer.csv')

pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pf.set_control(control)
states = pf.groupby_state()
series = states['SecB WT apo']

temperature, pH = 273.15 + 30, 8.
kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

if guess:
    wt_avg_result = kf.weighted_avg_fit()
    output = wt_avg_result.output
    output.to_file(directory / 'test_data' / 'ecSecB_guess.txt')
else:
    output = txt_to_protein(directory / 'test_data' / 'ecSecB_guess.txt')


fr_torch = kf.global_fit(output, epochs=epochs)
fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.txt')
with open(directory / 'test_data' / 'ecSecB_torch_fit.pick', 'wb') as f:
    pickle.dump(fr_torch, f)

series_dimer = states['SecB his dimer apo']
kf_dimer = KineticsFitting(series_dimer, bounds=(1e-2, 800), temperature=temperature, pH=pH)
bf = BatchFitting([kf, kf_dimer], [output, output])

batch_result = bf.global_fit(epochs=epochs)

batch_result.output.to_csv(directory / 'test_data' / 'ecSecB_batch.csv')