from pathlib import Path
import numpy as np
from pyhdx import PeptideMasterTable, read_dynamx, HDXMeasurement

current_dir = Path(__file__).parent
np.random.seed(43)

fpath = current_dir.parent / 'tests' / 'test_data' / 'ecSecB_apo.csv'
data = read_dynamx(fpath)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(('Full deuteration control', 0.167))

sequence = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'

hdxm = HDXMeasurement(pmt.get_state('SecB WT apo'), sequence=sequence)
print(hdxm.coverage.protein)

#hdxm.coverage.protein.to_file('test.txt', fmt='pprint')