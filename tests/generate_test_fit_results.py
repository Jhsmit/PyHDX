from pyhdx.fileIO import read_dynamx, csv_to_protein
from pyhdx import PeptideMasterTable, KineticsFitting, BatchFitting, KineticsSeries
from pyhdx.models import HDXMeasurementSet
from pyhdx.fitting import fit_rates_weighted_average, fit_gibbs_global, fit_gibbs_global_batch
from pyhdx.local_cluster import default_client
import numpy as np
import pickle
from pathlib import Path
import torch

"""Run this file to renew the fit results which is used to test against"""

torch.manual_seed(43)
np.random.seed(43)
epochs = 1000
sequence =       'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'
sequence_dimer = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'


directory = Path(__file__).parent
test_data_dir = directory / 'test_data'

guess = False  # guess true requires dask cluster at config defined ip/port
control = ('Full deuteration control', 0.167)

data = read_dynamx(test_data_dir / 'ecSecB_apo.csv', test_data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(control)
temperature, pH = 273.15 + 30, 8.

series = KineticsSeries(pmt.get_state('SecB WT apo'), sequence=sequence, temperature=temperature, pH=pH)
kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

if guess:
    client = default_client()
    wt_avg_result = fit_rates_weighted_average(series, bounds=(1e-2, 800)).compute()
    output = wt_avg_result.output
    output.to_file(directory / 'test_data' / 'ecSecB_guess.txt')
else:
    output = csv_to_protein(directory / 'test_data' / 'ecSecB_guess.txt')

gibbs_guess = series.guess_deltaG(output['rate'])
fr_torch = fit_gibbs_global(series, gibbs_guess, epochs=epochs)
fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.txt')

series_dimer = KineticsSeries(pmt.get_state('SecB his dimer apo'), sequence=sequence_dimer,
                              temperature=temperature, pH=pH)

hdx_set = HDXMeasurementSet([series_dimer, series])

gibbs_guess = hdx_set.guess_deltaG([output['rate'], output['rate']])
batch_result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=epochs)

batch_result.output.to_file(directory / 'test_data' / 'ecSecB_batch.csv')
batch_result.output.to_file(directory / 'test_data' / 'ecSecB_batch.txt', fmt='pprint')

# Order is inverted compared to test!
mock_alignment = {
    'dimer':   'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------',
    'apo':     'MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA',
}
Nr_a = len(next(iter(mock_alignment.values())))
alignment_array = np.empty((2, Nr_a), dtype='U1')  # Only max two aligned proteins can be fitted globally together
alignment_array[0] = list(mock_alignment['dimer'])
alignment_array[1] = list(mock_alignment['apo'])


kf1 = KineticsFitting(series_dimer, pH=8, temperature=273.15+30)
kf2 = KineticsFitting(series, pH=8, temperature=273.15+30)

bf = BatchFitting([kf1, kf2], [output, output])

aligned_result = bf.global_fit_aligned(alignment_array, r1=2, r2=5, epochs=1000)

aligned_result.output.to_file(directory / 'test_data' / 'ecSecB_batch_aligned.csv')
aligned_result.output.to_file(directory / 'test_data' / 'ecSecB_batch_aligned.txt', fmt='pprint')

# kf_dimer = KineticsFitting(series_dimer, bounds=(1e-2, 800), temperature=temperature, pH=pH)
#
# bf = BatchFitting([kf, kf_dimer], [output, output])
# batch_result = bf.global_fit(epochs=epochs)

series.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.csv')
series.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.txt', fmt='pprint')