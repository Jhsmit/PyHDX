"""Load a HDX-MS dataset (peptide list with D-uptake per peptide, csv format)"""

from pathlib import Path
import numpy as np
from pyhdx import PeptideMasterTable, read_dynamx, HDXMeasurement

current_dir = Path(__file__).parent
output_dir = current_dir / 'output'
output_dir.mkdir(exist_ok=True)
np.random.seed(43)

fpath = current_dir.parent / 'tests' / 'test_data' / 'ecSecB_apo.csv'
data = read_dynamx(fpath)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(('Full deuteration control', 0.167))

sequence = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'

temperature, pH = 273.15 + 30, 8.
hdxm = HDXMeasurement(pmt.get_state('SecB WT apo'), sequence=sequence, pH=pH, temperature=temperature)

hdxm.coverage.protein.to_file(output_dir / 'SecB_info.txt', fmt='pprint')