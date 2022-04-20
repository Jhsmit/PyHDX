"""Obtain ΔG for ecSecB tetramer and dimer"""
from pathlib import Path
from pyhdx.batch_processing import yaml_to_hdxmset
from pyhdx.fileIO import csv_to_dataframe, save_fitresult
from pyhdx.fitting import fit_gibbs_global_batch
import yaml

cwd = Path(__file__).parent

data_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"

yaml_dict = yaml.safe_load(Path(data_dir / "data_states.yaml").read_text())

hdx_set = yaml_to_hdxmset(yaml_dict, data_dir=data_dir)

initial_guess_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

guesses = hdx_set[0].guess_deltaG(initial_guess_rates["rate"])
fit_kwargs = yaml.safe_load(Path(data_dir / "fit_settings.yaml").read_text())

fr = fit_gibbs_global_batch(hdx_set, guesses, **fit_kwargs)
save_fitresult(output_dir / "ecsecb_tetramer_dimer", fr)
