"""Load SecB HDX-MS data and guesses, perform global fit of gibbs free energy, save results to directory"""

# %%
from pathlib import Path

import numpy as np
import proplot as pplt
import yaml

from pyhdx.datasets import HDXDataSet
from pyhdx.fileIO import csv_to_dataframe, save_fitresult
from pyhdx.fitting import fit_gibbs_global, fit_rates_weighted_average
from pyhdx.local_cluster import default_client
from pyhdx.models import HDXMeasurement

# %%

guess = False  # Set to True to redo initial guesses
fit_kwargs = {"epochs": 200000, "lr": 1e4, "stop_loss": 1e-6}

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

# %%
# evaluate the fit, plot measured and fitted D-uptake for a few peptides
NUM_EVAL_POINTS = 100
all_timepoints = np.concatenate([hdxm.timepoints for hdxm in fr_torch.hdxm_set])

elem = all_timepoints[np.nonzero(all_timepoints)]
start = np.log10(elem.min())
end = np.log10(elem.max())
pad = (end - start) * 0.1
time = np.logspace(start - pad, end + pad, num=NUM_EVAL_POINTS, endpoint=True)


# %%
# evaluate the fitted model at timepoints
d_calc = fr_torch(time)  # Ns x Np x Nt
d_calc
# %%
# choose which sample to plot
# here we plot the first one (= SecB tetramer)
sample_idx = 0
hdxm = fr_torch.hdxm_set.hdxm_list[sample_idx]
d_calc_s = d_calc[sample_idx]

# %%
# make a subplot grid to plot the first 24 peptides
fig, axes = pplt.subplots(ncols=4, nrows=6)
for i, ax in enumerate(axes):
    ax.plot(time, d_calc_s[i], color="r")
    ax.scatter(hdxm.timepoints, hdxm.d_exp.iloc[i], color="k")

axes.format(xscale="log", xlabel="Time (s)", ylabel="Corrected D-uptake", ylim=(0, None))
# %%
