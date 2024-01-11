"""Load two HDX-MS datasets and guesses and perform fitting in batch with a second regualizer

Example of using checkpoints to save model history

"""
# %%

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pyhdx.fileIO import (
    csv_to_protein,
    read_dynamx,
    dataframe_to_file,
    save_fitresult,
)
from pyhdx.fitting import fit_gibbs_global_batch
from pyhdx.fitting_torch import CheckPoint
from pyhdx.models import PeptideMasterTable, HDXMeasurement, HDXMeasurementSet


# %%

# Pycharm scientific mode
if "__file__" not in locals():
    __file__ = Path().cwd() / "templates" / "script.py"

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
data_dir = current_dir.parent / "tests" / "test_data"

# %%

data = read_dynamx(data_dir / "input" / "ecSecB_apo.csv", data_dir / "input" / "ecSecB_dimer.csv")
pmt = PeptideMasterTable(data)
pmt.set_control(("Full deuteration control", 0.167 * 60))

dimer = HDXMeasurement(pmt.get_state("SecB his dimer apo"), pH=8, temperature=273.15 + 30)
apo = HDXMeasurement(pmt.get_state("SecB WT apo"), pH=8, temperature=273.15 + 30)

hdx_set = HDXMeasurementSet([dimer, apo])
guess = csv_to_protein(data_dir / "output" / "ecSecB_guess.csv")

rates_df = pd.DataFrame({name: guess["rate"] for name in hdx_set.names})
gibbs_guess = hdx_set.guess_deltaG(rates_df)


# %%
# Example fit with only 5000 epochs and high learning rate
# Checkpoint stores model history every `epoch_step` epochs
checkpoint = CheckPoint(epoch_step=250)
result = fit_gibbs_global_batch(
    hdx_set, gibbs_guess, r1=0.5, r2=0.1, epochs=5000, lr=1e5, callbacks=[checkpoint]
)
print(
    f"MSE loss: {result.mse_loss:.2f}, "
    f"Reg loss: {result.reg_loss:.2f}, "
    f"Reg percent: {result.regularization_percentage:.0f}%"
)


df = checkpoint.to_dataframe(hdx_set.names)
dataframe_to_file(output_dir / "model_history.csv", df)
dataframe_to_file(output_dir / "model_history.txt", df, fmt="pprint")


# Checkpoint history scatter plot
# Note that these are raw dG values including interpolated values in regions of no coverage
history = checkpoint.model_history
num = len(history)
cmap = mpl.cm.get_cmap("winter")
norm = mpl.colors.Normalize(vmin=1, vmax=num * checkpoint.epoch_step)
colors = iter(cmap(np.linspace(0, 1, num=num)))

fig, ax = plt.subplots()
for key, val in history.items():
    n = len(val["dG"].numpy().squeeze())
    ax.scatter(hdx_set.coverage.index, val["dG"].numpy().squeeze()[0], color=next(colors))

fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm), label="Epochs")
plt.show()

# Human readable output
result.to_file(output_dir / "Batch_fit_result.txt", fmt="pprint")

# Machine readable output
result.to_file(output_dir / "Batch_fit_result.csv", fmt="csv")

# Save full fitresult
save_fitresult(output_dir / "SecB_tetramer_dimer_batch", result)
