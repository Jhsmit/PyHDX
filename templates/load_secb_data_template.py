from pathlib import Path
import numpy as np
from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries

current_dir = Path(__file__).parent
np.random.seed(43)

fpath = current_dir.parent / 'tests' / 'test_data' / 'ecSecB_apo.csv'
data = read_dynamx(fpath)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(('Full deuteration control', 0.167))

sequence = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'

series = KineticsSeries(pmt.get_state('SecB WT apo'), sequence=sequence)
print(series)

#series.coverage.protein.to_file('test.txt', fmt='pprint')

from pyhdx.fileIO import csv_to_protein
protein = csv_to_protein(current_dir.parent / 'tests' / 'test_data' / 'ecSecB_info.csv', column_depth=1)

print(protein.df)
#print(protein.index)