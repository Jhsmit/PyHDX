"""Load SecB HDX-MS data and guesses, perform global fit of gibbs free energy, save results to directory"""

# %%
from pathlib import Path

from pyhdx.datasets import HDXDataSet
from pyhdx.models import HDXMeasurement
from pyhdx.fitting import fit_gibbs_global, fit_rates_weighted_average
from pyhdx.fileIO import csv_to_dataframe, save_fitresult
from pyhdx.local_cluster import default_client
import yaml

# %%

guess = False  # Set to True to redo initial guesses
fit_kwargs = {"epochs": 10000, "lr": 1e4, "stop_loss": 1e-6}

# %%
current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
hdx_spec = yaml.safe_load(yaml_stream)

# %%

input_dir = current_dir.parent / "tests" / "test_data" / "input"
dataset = HDXDataSet.from_spec(hdx_spec, data_dir=input_dir)

# %%
hdxm = HDXMeasurement.from_dataset(dataset, state="SecB_tetramer")
print(hdxm.timepoints)

# %%

if guess:
    client = default_client()
    wt_avg_result = fit_rates_weighted_average(hdxm, client=client)
    init_guess = wt_avg_result.output
else:
    init_guess = csv_to_dataframe(
        current_dir.parent / "tests" / "test_data" / "output" / "ecSecB_guess.csv"
    )

gibbs_guess = hdxm.guess_deltaG(init_guess["rate"])

# %%

fr_torch = fit_gibbs_global(hdxm, gibbs_guess, **fit_kwargs)

# Human readable output
fr_torch.to_file(output_dir / "SecB_fit_result.txt", fmt="pprint")

# Machine readable output
fr_torch.to_file(output_dir / "SecB_fit_result.csv", fmt="csv")

save_fitresult(output_dir / "SecB_fit", fr_torch)
