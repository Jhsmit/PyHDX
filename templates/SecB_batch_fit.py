from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, KineticsSeries
from pyhdx.fitting import BatchFitting, KineticsFitting
from pyhdx.fileIO import csv_to_protein

current_dir = Path(__file__).parent

data_dir = current_dir.parent / 'tests' / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))

st1 = KineticsSeries(pmt.get_state('SecB his dimer apo'))
st2 = KineticsSeries(pmt.get_state('SecB WT apo'))

kf1 = KineticsFitting(st1, pH=8, temperature=273.15+30)
kf2 = KineticsFitting(st2, pH=8, temperature=273.15+30)

guesses = csv_to_protein(data_dir / 'ecSecB_guess.txt')

bf = BatchFitting([kf1, kf2], [guesses, guesses])

result = bf.global_fit(epochs=10)

#Human readable output
result.output.to_file('Batch_fit_result.txt', fmt='pprint')

#Machine readable output
result.output.to_file('Batch_fit_result.csv', fmt='csv')