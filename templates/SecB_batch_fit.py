from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, KineticsSeries, HDXMeasurementSet
from pyhdx.fitting import BatchFitting, fit_gibbs_global_batch
from pyhdx.fileIO import csv_to_protein

current_dir = Path(__file__).parent

data_dir = current_dir.parent / 'tests' / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))

st1 = KineticsSeries(pmt.get_state('SecB his dimer apo'), pH=8, temperature=273.15+30)
st2 = KineticsSeries(pmt.get_state('SecB WT apo'), pH=8, temperature=273.15+30)

hdx_set = HDXMeasurementSet([st1, st2])
guess = csv_to_protein(data_dir / 'ecSecB_guess.txt')

gibbs_guess = hdx_set.guess_deltaG([guess['rate'], guess['rate']])
result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=1000)

#Human readable output
result.output.to_file('Batch_fit_result.txt', fmt='pprint')

#Machine readable output
result.output.to_file('Batch_fit_result.csv', fmt='csv')