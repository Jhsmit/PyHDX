"""
Automagically plot all available figures from a fit result
"""

from pyhdx.fileIO import load_fitresult
from pyhdx.plot import linear_bars_figure
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# %%

__file__ = Path().cwd() / "templates" / "script.py"  # Uncomment for PyCharm scientific mode

cwd = Path(__file__).parent
output_dir = cwd / "output" / "single_figures"
output_dir.mkdir(exist_ok=True)
fit_result = load_fitresult(cwd / "output" / "SecB_tetramer_dimer_batch")

# rfus = fit_result.hdxm_set.rfu_residues
# %%
# self = fit_result.hdxm_set
# [hdxm.rfu_residues for hdxm in self]
# rfus = pd.concat([hdxm.rfu_residues for hdxm in self],
#           keys=self.names, names=['state', 'exposure'], axis=1)
# add quantity multiindex level

rfus = fit_result.hdxm_set.rfu_residues
columns = pd.MultiIndex.from_tuples(
    [(*tup, "rfu") for tup in rfus.columns], names=[*rfus.columns.names, "quantity"]
)
rfus.columns = columns


# %%

# RFUs grouped by state
name = "linear_bars_rfu_by_state"
fig, axes, cbar = linear_bars_figure(rfus, field="rfu")
plt.savefig(output_dir / f"{name}.png")
plt.savefig(output_dir / f"{name}.pdf")

# %%

# RFUs grouped by exposure
name = "linear_bars_rfu_by_exposure"
fig, axes, cbar = linear_bars_figure(rfus, field="rfu", groupby="exposure")
plt.savefig(output_dir / f"{name}.png")
plt.savefig(output_dir / f"{name}.pdf")

# %%

# dRFUs grouped by exposure, wrt first sample (test - reference, positive dRFU is more flexible)
name = "linear_bars_drfu_by_exposure"
fig, axes, cbar = linear_bars_figure(rfus, field="rfu", groupby="exposure", reference=0)
plt.savefig(output_dir / f"{name}.png")
plt.savefig(output_dir / f"{name}.pdf")


# %%

# dGs
name = "linear_bars_dG"
fig, axes, cbar = linear_bars_figure(fit_result.output, field="dG")
plt.savefig(output_dir / f"{name}.png")
plt.savefig(output_dir / f"{name}.pdf")

# %%

# ddGs
name = "linear_bars_ddG"
fig, axes, cbar = linear_bars_figure(fit_result.output, field="dG", reference=0)
plt.savefig(output_dir / f"{name}.png")
plt.savefig(output_dir / f"{name}.pdf")
