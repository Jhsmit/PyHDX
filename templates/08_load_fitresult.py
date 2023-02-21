"""Load a fit result from a directory where a fit result was saved with save_fitresult"""
from pyhdx.fileIO import load_fitresult
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# %%

current_dir = Path(__file__).parent
print(current_dir.parent / "tests" / "test_data" / "output" / "ecsecb_tetramer_dimer")

fit_result = load_fitresult(
    current_dir.parent / "tests" / "test_data" / "output" / "ecsecb_tetramer_dimer"
)

# %%

mse_df = fit_result.get_peptide_mse()


s = 0  # index of protein state to view
p = 20  # index of the peptide to view

timepoints = fit_result.hdxm_set[s].timepoints
timepoints = timepoints[np.nonzero(timepoints)]
tmin = np.log10(timepoints.min())
tmax = np.log10(timepoints.max())

pad = 0.1 * (tmax - tmin)

time = np.logspace(tmin - pad, tmax + pad, num=100)

d_calc = fit_result(time)
d_exp = fit_result.hdxm_set.d_exp

fit_result.losses[["mse_loss", "reg_1"]].plot()

fig, ax = plt.subplots()
ax.scatter(fit_result.hdxm_set[s].timepoints, d_exp[s, p, :], color="k")
ax.plot(time, d_calc[s, p, :], color="r")
ax.set_xscale("log")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Corrected D-uptake (Da)")
ax.set_title(f"Peptide index: {p}")

output = fit_result.output
protein_name = output.columns.unique(level=0)[0]

fig, ax = plt.subplots()
ax.set_aspect(2)
ax.scatter(fit_result.output.index, fit_result.output[protein_name]["dG"] * 1e-3)
ax.set_xlabel("Residue number")
ax.set_ylabel("ΔG (kJ/mol)")
plt.show()
