from pathlib import Path

import numpy as np
import pandas as pd
import torch
import yaml

from pyhdx import HDXMeasurement
from pyhdx.batch_processing import StateParser
from pyhdx.fileIO import (
    read_dynamx,
    save_fitresult,
    dataframe_to_file,
    csv_to_dataframe,
)
from pyhdx.fitting import (
    fit_rates_weighted_average,
    fit_gibbs_global,
    fit_gibbs_global_batch,
    fit_gibbs_global_batch_aligned,
    fit_d_uptake,
)
from pyhdx.models import HDXMeasurementSet
from pyhdx.process import apply_control, correct_d_uptake, filter_peptides

"""Run this file to renew the fit results which is used to test against"""

# Toggle to also generate long computation time fits
do_long_fit = False
epochs_long = 20000

torch.manual_seed(43)
np.random.seed(43)
epochs = 1000
sequence = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"
sequence_dimer = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"

guess = False

df_apo = read_dynamx(input_dir / "ecSecB_apo.csv")
df_dimer = read_dynamx(input_dir / "ecSecB_dimer.csv")

fd = {"state": "Full deuteration control", "exposure": {"value": 0.167, "unit": "min"}}
fd_df = filter_peptides(df_apo, **fd)

apo_peptides = filter_peptides(df_apo, state="SecB WT apo")
dimer_peptides = filter_peptides(df_dimer, state="SecB his dimer apo")

apo_control = apply_control(apo_peptides, fd_df)
dimer_control = apply_control(dimer_peptides, fd_df)

apo_corrected = correct_d_uptake(apo_control)
dimer_corrected = correct_d_uptake(dimer_control)

temperature, pH = 273.15 + 30, 8.0

hdxm_apo = HDXMeasurement(
    apo_corrected, sequence=sequence, temperature=temperature, pH=pH
)

reduced_data = apo_corrected[apo_corrected["end"] < 40]
hdxm_reduced = HDXMeasurement(reduced_data, temperature=temperature, pH=pH, c_term=155)

result = fit_rates_weighted_average(hdxm_reduced)
reduced_guess = result.output
dataframe_to_file(output_dir / "ecSecB_reduced_guess.csv", reduced_guess)
dataframe_to_file(output_dir / "ecSecB_reduced_guess.txt", reduced_guess, fmt="pprint")

gibbs_guess = hdxm_reduced.guess_deltaG(reduced_guess["rate"])
fr_torch = fit_gibbs_global(hdxm_reduced, gibbs_guess, epochs=epochs, r1=2)
save_fitresult(output_dir / "ecsecb_reduced", fr_torch)

if guess:
    wt_avg_result = fit_rates_weighted_average(
        hdxm_apo, bounds=(1e-2 / 60.0, 800 / 60.0)
    )
    guess_output = wt_avg_result.output
    dataframe_to_file(output_dir / "ecSecB_guess.csv", guess_output)
    dataframe_to_file(output_dir / "ecSecB_guess.txt", guess_output, fmt="pprint")
else:
    guess_output = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

# Export protein sequence and intrinsic rate of exchange
hdxm_apo.coverage.protein.to_file(output_dir / "ecSecB_info.csv")
hdxm_apo.coverage.protein.to_file(output_dir / "ecSecB_info.txt", fmt="pprint")

# Save RFU values
rfu_df = hdxm_apo.rfu_residues
dataframe_to_file(output_dir / "ecSecB_rfu_per_exposure.csv", rfu_df)
dataframe_to_file(output_dir / "ecSecB_rfu_per_exposure.txt", rfu_df, fmt="pprint")

# Save data DataFrame
data_df = hdxm_apo.data
dataframe_to_file(output_dir / "ecSecB_data.csv", data_df)
dataframe_to_file(output_dir / "ecSecB_data.txt", data_df, fmt="pprint")

# Fit D-uptake per timepoint
fr_d = fit_d_uptake(hdxm_apo, r1=1.0, repeats=3)
dataframe_to_file(output_dir / f"ecSecB_d_uptake.csv", fr_d.output)

gibbs_guess = hdxm_apo.guess_deltaG(guess_output["rate"])
fr_torch = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=epochs, r1=2)

dataframe_to_file(output_dir / f"ecSecB_torch_fit.csv", fr_torch.output)
dataframe_to_file(output_dir / f"ecSecB_torch_fit.txt", fr_torch.output, fmt="pprint")

fr_torch = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=epochs_long, r1=2)
dataframe_to_file(
    output_dir / f"ecSecB_torch_fit_epochs_{epochs_long}.csv", fr_torch.output
)
dataframe_to_file(
    output_dir / f"ecSecB_torch_fit_epochs_{epochs_long}.txt",
    fr_torch.output,
    fmt="pprint",
)

# ----------
# Batch fits
# ----------

hdxm_dimer = HDXMeasurement(
    dimer_corrected,
    sequence=sequence_dimer,
    temperature=temperature,
    pH=pH,
)

hdx_set = HDXMeasurementSet([hdxm_dimer, hdxm_apo])

rates_df = pd.DataFrame({name: guess_output["rate"] for name in hdx_set.names})
gibbs_guess = hdx_set.guess_deltaG(rates_df)
batch_result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=epochs)

dataframe_to_file(output_dir / "ecSecB_batch.csv", batch_result.output)
dataframe_to_file(output_dir / "ecSecB_batch.txt", batch_result.output, fmt="pprint")

# Save errors and losses
dataframe_to_file(
    output_dir / "ecSecB_batch_peptide_mse.csv", batch_result.get_peptide_mse()
)
dataframe_to_file(
    output_dir / "ecSecB_batch_residue_mse.csv", batch_result.get_residue_mse()
)
dataframe_to_file(output_dir / "ecSecB_batch_loss.csv", batch_result.losses)

# Aligned sequences test
# -------------

# Order is inverted compared to test!
mock_alignment = {
    "dimer": "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------",
    "apo": "MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA",
}

hdx_set.add_alignment(list(mock_alignment.values()))

aligned_result = fit_gibbs_global_batch_aligned(
    hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000
)

dataframe_to_file(output_dir / "ecSecB_batch_aligned.csv", aligned_result.output)
dataframe_to_file(
    output_dir / "ecSecB_batch_aligned.txt", aligned_result.output, fmt="pprint"
)

# ------------------
# Reduced Batch fits
# ------------------

dimer_red = dimer_corrected[dimer_corrected["end"] < 40]
hdxm_reduced_dimer = HDXMeasurement(
    dimer_red,
    temperature=temperature,
    pH=pH,
    c_term=155,
)

reduced_hdx_set = HDXMeasurementSet([hdxm_reduced, hdxm_reduced_dimer])
gibbs_guess = reduced_hdx_set[0].guess_deltaG(reduced_guess["rate"])
batch_result = fit_gibbs_global_batch(reduced_hdx_set, gibbs_guess, epochs=epochs)
save_fitresult(output_dir / "ecsecb_reduced_batch", batch_result)

# -------------------
# delta N/C tail fits
# -------------------

# These datasets have unequal overall coverage range between datasets
yaml_file = input_dir / "data_states_deltas.yaml"
yaml_dict = yaml.safe_load(yaml_file.read_text())
parser = StateParser(yaml_dict, data_src=input_dir)
hdxm_set = parser.load_hdxmset()
# hdxm_set = yaml_to_hdxmset(yaml_dict, data_dir=input_dir)

gibbs_guess = hdxm_set[0].guess_deltaG(guess_output["rate"])

batch_result = fit_gibbs_global_batch(hdxm_set, gibbs_guess, epochs=200)
save_fitresult(output_dir / "ecsecb_delta_batch", batch_result)
