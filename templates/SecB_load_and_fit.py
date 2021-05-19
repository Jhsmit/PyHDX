from pathlib import Path
import numpy as np
from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries
from pyhdx.fitting import fit_gibbs_global, fit_rates_weighted_average
from pyhdx.fileIO import csv_to_protein
from pyhdx.local_cluster import default_client
import asyncio

guess = False
epochs = 1000
root_dir = Path().resolve().parent
test_data_dir = root_dir / 'tests' / 'test_data'
input_file_path = test_data_dir / 'ecSecB_apo.csv'

# Load the data of two Dynamx files, and combine the result to one table
data = read_dynamx(test_data_dir / 'ecSecB_apo.csv', test_data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(('Full deuteration control', 0.167))
temperature, pH = 273.15 + 30, 8.
series = KineticsSeries(pmt.get_state('SecB WT apo'), temperature=temperature, pH=pH)


#kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

if guess:
    client = default_client()
    wt_avg_result = fit_rates_weighted_average(series).compute()
    init_guess = wt_avg_result.output
else:
    init_guess = csv_to_protein(test_data_dir / 'ecSecB_guess.txt')

gibbs_guess = series.guess_deltaG(init_guess['rate'])
fr_torch = fit_gibbs_global(series, gibbs_guess, epochs=epochs)
print(fr_torch.metadata['total_loss'])

