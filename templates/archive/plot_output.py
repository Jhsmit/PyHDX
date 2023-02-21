"""
Automagically plot all available figures from a fit result
"""

from pyhdx.fileIO import load_fitresult
from pyhdx.plot import FitResultPlot
from pathlib import Path

# %%

# __file__ = Path().cwd() / 'templates'/ 'script.py'  # Uncomment for PyCharm scientific mode


cwd = Path(__file__).parent
output_dir = cwd / "output" / "figures"
output_dir.mkdir(exist_ok=True)
fit_result = load_fitresult(cwd / "output" / "SecB_tetramer_dimer_batch")

fr_plot = FitResultPlot(fit_result, output_path=output_dir)

kwargs = {
    "residue_scatter": {"cmap": "BuGn"},  # change default colormap
    "ddG_scatter": {
        "reference": 1
    },  # Set reference for ΔΔG to the second (index 1 state) (+ APO state (tetramer))
}

fr_plot.plot_all(**kwargs)
