from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable
from pyhdx.fitting import BatchFitting, KineticsFitting
from pyhdx.fileIO import txt_to_protein, csv_to_protein
import asyncio

current_dir = Path(__file__).parent

data_dir = current_dir.parent / 'tests' / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))
states = pmt.groupby_state()

st1 = states['SecB his dimer apo']
st2 = states['SecB WT apo']

kf1 = KineticsFitting(st1, pH=8, temperature=273.15+30)
kf2 = KineticsFitting(st2, pH=8, temperature=273.15+30)

print(kf1.cluster)
if __name__ == '__main__':
    guesses = csv_to_protein(data_dir / 'ecSecB_guess.txt')
    result = asyncio.run(kf1.global_fit_async(guesses, epochs=10))
    print(result.output)

