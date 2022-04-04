from pyhdx.fileIO import read_dynamx, save_fitresult, dataframe_to_file, csv_to_dataframe
from pyhdx import PeptideMasterTable, HDXMeasurement
from pyhdx.models import HDXMeasurementSet
from pyhdx.fitting import fit_rates_weighted_average, fit_gibbs_global, fit_gibbs_global_batch, fit_gibbs_global_batch_aligned
from pyhdx.local_cluster import default_client
from pyhdx.batch_processing import yaml_to_hdxmset
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
import torch


torch.manual_seed(43)
np.random.seed(43)
epochs = 1000
sequence =       'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'
sequence_dimer = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'


cwd = Path(__file__).parent.parent
input_dir = cwd / 'test_data' / 'input'
output_dir = cwd / 'test_data' / 'output'

control = ('Full deuteration control', 0.167*60)

data = read_dynamx(input_dir / 'ecSecB_apo.csv', input_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(control)
temperature, pH = 273.15 + 30, 8.

hdxm_apo = HDXMeasurement(pmt.get_state('SecB WT apo'), temperature=temperature, pH=pH)
hdxm_dimer = HDXMeasurement(pmt.get_state('SecB his dimer apo'),
                            temperature=temperature, pH=pH)

guess = csv_to_dataframe(output_dir / 'ecSecB_guess.csv')


gibbs_guess = hdxm_apo.guess_deltaG(guess['rate'])
fr_torch = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=epochs, r1=2)

# ----------
# Batch fits
# ----------


hdx_set = HDXMeasurementSet([hdxm_dimer, hdxm_apo])

rates_df = pd.DataFrame({name: guess['rate'] for name in hdx_set.names})
gibbs_guess = hdx_set.guess_deltaG(rates_df)
batch_result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=epochs)

dataframe_to_file(output_dir / 'ecSecB_batch.csv', batch_result.output)
dataframe_to_file(output_dir / 'ecSecB_batch.txt', batch_result.output, fmt='pprint')

# Save errors and losses
dataframe_to_file(output_dir / 'ecSecB_batch_peptide_mse.csv', batch_result.get_peptide_mse())
dataframe_to_file(output_dir / 'ecSecB_batch_residue_mse.csv', batch_result.get_residue_mse())
dataframe_to_file(output_dir / 'ecSecB_batch_loss.csv', batch_result.losses)

