from pyhdx.fileIO import read_dynamx, csv_to_protein
from pyhdx import PeptideMasterTable, HDXMeasurement
from pyhdx.models import HDXMeasurementSet
from pyhdx.fitting import fit_rates_weighted_average, fit_gibbs_global, fit_gibbs_global_batch, fit_gibbs_global_batch_aligned
from pyhdx.local_cluster import default_client
import numpy as np
from pathlib import Path
import torch

"""Run this file to renew the fit results which is used to test against"""

# Toggle to also generate long computation time fits
do_long_fit = False
epochs_long = 20000

torch.manual_seed(43)
np.random.seed(43)
epochs = 1000
sequence =       'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'
sequence_dimer = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'


directory = Path(__file__).parent
test_data_dir = directory / 'test_data'

guess = True  # guess true requires dask cluster at config defined ip/port
control = ('Full deuteration control', 0.167)

data = read_dynamx(test_data_dir / 'ecSecB_apo.csv', test_data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(control)
temperature, pH = 273.15 + 30, 8.

hdxm = HDXMeasurement(pmt.get_state('SecB WT apo'), sequence=sequence, temperature=temperature, pH=pH)

data = pmt.get_state('SecB WT apo')
reduced_data = data[data['end'] < 40]
reduced_hdxm = HDXMeasurement(reduced_data)

result = fit_rates_weighted_average(reduced_hdxm)
output = result.output
output.to_file(directory / 'test_data' / 'ecSecB_reduced_guess.csv')
output.to_file(directory / 'test_data' / 'ecSecB_reduced_guess.txt', fmt='pprint')

if guess:
    # Initial guesses changed (by refactoring to rfu?)
    from dask.distributed import Client
    client = Client('tcp://127.0.0.1:52348')
    wt_avg_result = fit_rates_weighted_average(hdxm, bounds=(1e-2, 800))
    output = wt_avg_result.output
    output.to_file(directory / 'test_data' / 'ecSecB_guess.csv')
else:
    output = csv_to_protein(directory / 'test_data' / 'ecSecB_guess.csv')
    # import pandas as pd
    # output = pd.read_csv(directory / 'test_data' / 'ecSecB_guess.txt', header=[0], comment='#', index_col=0)

gibbs_guess = hdxm.guess_deltaG(output['rate'])
fr_torch = fit_gibbs_global(hdxm, gibbs_guess, epochs=epochs, r1=2)
fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.csv')
fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.txt', fmt='pprint')

fr_torch = fit_gibbs_global(hdxm, gibbs_guess, epochs=epochs_long, r1=2)
fr_torch.output.to_file(directory / 'test_data' / f'ecSecB_torch_fit_epochs_{epochs_long}.csv')
fr_torch.output.to_file(directory / 'test_data' / f'ecSecB_torch_fit_epochs_{epochs_long}.txt', fmt='pprint')


# ---------------
# Batch fits

hdxm_dimer = HDXMeasurement(pmt.get_state('SecB his dimer apo'), sequence=sequence_dimer,
                            temperature=temperature, pH=pH)

hdx_set = HDXMeasurementSet([hdxm_dimer, hdxm])

gibbs_guess = hdx_set.guess_deltaG([output['rate'], output['rate']])
batch_result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=epochs)

batch_result.output.to_file(directory / 'test_data' / 'ecSecB_batch.csv')
batch_result.output.to_file(directory / 'test_data' / 'ecSecB_batch.txt', fmt='pprint')

# Order is inverted compared to test!
mock_alignment = {
    'dimer':   'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------',
    'apo':     'MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA',
}

hdx_set.add_alignment(list(mock_alignment.values()))

aligned_result = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000)

aligned_result.output.to_file(directory / 'test_data' / 'ecSecB_batch_aligned.csv')
aligned_result.output.to_file(directory / 'test_data' / 'ecSecB_batch_aligned.txt', fmt='pprint')

hdxm.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.csv')
hdxm.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.txt', fmt='pprint')