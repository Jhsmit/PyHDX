from pyhdx.fileIO import read_dynamx, csv_to_protein, save_fitresult
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

guess = True
control = ('Full deuteration control', 0.167*60)

data = read_dynamx(test_data_dir / 'ecSecB_apo.csv', test_data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(control)
temperature, pH = 273.15 + 30, 8.

hdxm = HDXMeasurement(pmt.get_state('SecB WT apo'), sequence=sequence, temperature=temperature, pH=pH)

data = pmt.get_state('SecB WT apo')
reduced_data = data[data['end'] < 40]
hdxm_reduced = HDXMeasurement(reduced_data, temperature=temperature, pH=pH)

result = fit_rates_weighted_average(hdxm_reduced)
reduced_guess = result.output
reduced_guess.to_file(directory / 'test_data' / 'ecSecB_reduced_guess.csv')
reduced_guess.to_file(directory / 'test_data' / 'ecSecB_reduced_guess.txt', fmt='pprint')

gibbs_guess = hdxm_reduced.guess_deltaG(reduced_guess['rate'])
print(gibbs_guess)
fr_torch = fit_gibbs_global(hdxm_reduced, gibbs_guess, epochs=epochs, r1=2)
save_fitresult(directory / 'test_data' / 'ecsecb_reduced', fr_torch)

if guess:
    wt_avg_result = fit_rates_weighted_average(hdxm, bounds=(1e-2/60., 800/60.))
    output = wt_avg_result.output
    output.to_file(directory / 'test_data' / 'ecSecB_guess.csv')
    output.to_file(directory / 'test_data' / 'ecSecB_guess.txt', fmt='pprint')

else:
    output = csv_to_protein(directory / 'test_data' / 'ecSecB_guess.csv')

hdxm.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.csv')
hdxm.coverage.protein.to_file(directory / 'test_data' / 'ecSecB_info.txt', fmt='pprint')

gibbs_guess = hdxm.guess_deltaG(output['rate'])
print(gibbs_guess)
fr_torch = fit_gibbs_global(hdxm, gibbs_guess, epochs=epochs, r1=2)
fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.csv')
fr_torch.output.to_file(directory / 'test_data' / 'ecSecB_torch_fit.txt', fmt='pprint')

fr_torch = fit_gibbs_global(hdxm, gibbs_guess, epochs=epochs_long, r1=2)
fr_torch.output.to_file(directory / 'test_data' / f'ecSecB_torch_fit_epochs_{epochs_long}.csv')
fr_torch.output.to_file(directory / 'test_data' / f'ecSecB_torch_fit_epochs_{epochs_long}.txt', fmt='pprint')


# ----------
# Batch fits
# ----------

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



# ------------------
# Reduced Batch fits
# ------------------

data = pmt.get_state('SecB his dimer apo')
reduced_data = data[data['end'] < 40]
hdxm_reduced_dimer = HDXMeasurement(pmt.get_state('SecB his dimer apo'), temperature=temperature, pH=pH)

reduced_hdx_set = HDXMeasurementSet([hdxm_reduced, hdxm_reduced_dimer])
gibbs_guess = reduced_hdx_set.guess_deltaG([reduced_guess['rate'], reduced_guess['rate']])
batch_result = fit_gibbs_global_batch(reduced_hdx_set, gibbs_guess, epochs=epochs)
save_fitresult(directory / 'test_data' / 'ecsecb_reduced_batch', batch_result)