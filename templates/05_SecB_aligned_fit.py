"""Load two HDX-MS datasets and guesses and perform fitting in batch with a second regualizer and mock alignment"""

from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, HDXMeasurement, HDXMeasurementSet
from pyhdx.fitting import fit_gibbs_global_batch_aligned
from pyhdx.fileIO import csv_to_protein
from pyhdx.alignment import align_dataframes
import numpy as np

mock_alignment = {
    'dimer':   'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------',
    'apo':     'MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA',
}

current_dir = Path(__file__).parent
output_dir = current_dir / 'output'
output_dir.mkdir(exist_ok=True)
data_dir = current_dir.parent / 'tests' / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))

st1 = HDXMeasurement(pmt.get_state('SecB his dimer apo'), pH=8, temperature=273.15 + 30)
st2 = HDXMeasurement(pmt.get_state('SecB WT apo'), pH=8, temperature=273.15 + 30)

guess = csv_to_protein(data_dir / 'ecSecB_guess.txt', header=[2], index_col=0)

hdx_set = HDXMeasurementSet([st1, st2])
gibbs_guess = hdx_set.guess_deltaG([guess['rate'], guess['rate']])
hdx_set.add_alignment(list(mock_alignment.values()))
result = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000)

#Human readable output
result.to_file(output_dir / 'Batch_aligned_fit_result.txt', fmt='pprint')

#Machine readable output
result.to_file(output_dir / 'Batch_aligned_fit_result.csv', fmt='csv')


