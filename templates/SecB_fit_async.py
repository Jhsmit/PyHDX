from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, KineticsSeries
from pyhdx.fitting import BatchFitting, KineticsFitting
from pyhdx.fileIO import csv_to_protein
import asyncio

# Requires Dask cluster at cfg adress to run

current_dir = Path(__file__).parent

data_dir = current_dir.parent / 'tests' / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))
#todo this should return dict of arrays or nested array


state_data = pmt.get_state('SecB his dimer apo')
state = KineticsSeries(state_data)

kf = KineticsFitting(state, pH=8, temperature=273.15+30)

if __name__ == '__main__':
    guesses = csv_to_protein(data_dir / 'ecSecB_guess.txt')
    result = asyncio.run(kf.global_fit_async(guesses, epochs=10))
    print(result.output)

