"""Load two HDX-MS datasets and guesses and perform fitting in batch with a secon regualizer"""
from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, HDXMeasurement, HDXMeasurementSet
from pyhdx.fitting import fit_gibbs_global_batch
from pyhdx.fileIO import csv_to_protein

current_dir = Path(__file__).parent
output_dir = current_dir / 'output'
output_dir.mkdir(exist_ok=True)
data_dir = current_dir.parent / 'tests' / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))

st1 = HDXMeasurement(pmt.get_state('SecB his dimer apo'), pH=8, temperature=273.15 + 30)
st2 = HDXMeasurement(pmt.get_state('SecB WT apo'), pH=8, temperature=273.15 + 30)

hdx_set = HDXMeasurementSet([st1, st2])
#todo update to new format
guess = csv_to_protein(data_dir / 'ecSecB_guess.txt', header=[2], index_col=0)

gibbs_guess = hdx_set.guess_deltaG([guess['rate'], guess['rate']])

# Example fit with only 1000 epochs and high regularizers
# For real data start with parameters r1=0.05, r2=0.5, epochs=100000
result = fit_gibbs_global_batch(hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000)
print(f"MSE loss: {result.mse_loss:.2f}, "
      f"Reg loss: {result.reg_loss:.2f}, "
      f"Reg percent: {result.regularization_percentage:.0f}%")

#Human readable output
result.to_file(output_dir / 'Batch_fit_result.txt', fmt='pprint')

#Machine readable output
result.to_file(output_dir / 'Batch_fit_result.csv', fmt='csv')