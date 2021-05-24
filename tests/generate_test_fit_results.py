from pyhdx.fileIO import read_dynamx, csv_to_protein
from pyhdx import PeptideMasterTable, KineticsFitting, BatchFitting, KineticsSeries
import numpy as np
import pickle
from pathlib import Path
import torch

"""Run this file to renew the fit results which is used to test against"""

torch.manual_seed(43)
np.random.seed(43)
epochs = 1000
sequence =       'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'
sequence_dimer = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'

print(sequence)

directory = Path(__file__).parent
test_data_dir = directory / 'test_data'

guess = False
control = ('Full deuteration control', 0.167)

data = read_dynamx(test_data_dir / 'ecSecB_apo.csv', test_data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(control)
series = KineticsSeries(pmt.get_state('SecB WT apo'), sequence=sequence)

temperature, pH = 273.15 + 30, 8.
kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

if guess:
    wt_avg_result = kf.weighted_avg_fit()
    output = wt_avg_result.output
    output.to_file(directory / 'test_data' / 'ecSecB_guess.txt')
else:
    output = csv_to_protein(directory / 'test_data' / 'ecSecB_guess.txt')

fr_torch = kf.global_fit(output, epochs=epochs)
temp = fr_torch.output

fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.txt')

series_dimer = KineticsSeries(pmt.get_state('SecB his dimer apo'), sequence=sequence_dimer)
kf_dimer = KineticsFitting(series_dimer, bounds=(1e-2, 800), temperature=temperature, pH=pH)
bf = BatchFitting([kf, kf_dimer], [output, output])

batch_result = bf.global_fit(epochs=epochs)

batch_result.output.to_file(directory / 'test_data' / 'ecSecB_batch.csv')
batch_result.output.to_file(directory / 'test_data' / 'ecSecB_batch.txt', fmt='pprint')

series.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.csv')
series.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.txt', fmt='pprint')