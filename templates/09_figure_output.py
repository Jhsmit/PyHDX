# %%
from pathlib import Path

import pandas as pd
import ultraplot as uplt

# from pymol import cmd
from pyhdx.config import cfg
from pyhdx.fileIO import csv_to_dataframe
from pyhdx.plot import (
    CMAP_NORM_DEFAULTS,
    ddG_scatter_figure,
    dG_scatter_figure,
    linear_bars_figure,
)
from pyhdx.support import apply_cmap, color_pymol

# %%
cwd = Path(__file__).parent
root_dir = cwd.parent
web_data_dir = root_dir / "tests" / "test_data" / "output" / "web"


# %%
dG_df = csv_to_dataframe(web_data_dir / "dG.csv")
plot_data = dG_df.xs("Gibbs_fit_1", axis=1, level=0)

# %%

width = cfg.plotting.page_width
aspect = cfg.plotting.dG_aspect
figure_kwargs = {
    "width": width,
    "refaspect": aspect,
}

figure_kwargs["ncols"] = 2

# %%
protein_states = plot_data.columns.get_level_values(0).unique()
protein_states

# %%

fig, axes, cbars = dG_scatter_figure(plot_data)  # , **figure_kwargs)
uplt.show()

# %%
# %%
linear_bars_figure(dG_df, groupby="fit_ID")
uplt.show()

# %%

linear_bars_figure(dG_df, groupby="fit_ID", reference="SecB_tetramer")
uplt.show()

# %%
protein_states = plot_data.columns.get_level_values(0).unique()
ddG_scatter_figure(plot_data, reference=protein_states[0])

# %%
from pymol import cmd

# Creating a colored structure
cmd.load("https://files.rcsb.org/download/1QYN.pdb")
cmd.set("antialias", 2)
cmd.set("fog", 0)

# %%
cmd.remove("resn HOH")  # This removes only water molecules
cmd.ipython_image()
# %%
# take dG values for SecB tetramer and reindex (pad with nan)
# such that these regions are correctly colored as no coverage
dG_values = plot_data[("SecB_tetramer", "dG")].reindex(pd.RangeIndex(1, 160))

# %%
# create a pandas Series with hexadeicmal codes
cmap, norm = CMAP_NORM_DEFAULTS["dG"]
colors = apply_cmap(dG_values, cmap, norm)
colors

# %%
# apply the colors to the pymol structure
color_pymol(colors, cmd)
cmd.ipython_image()

# %%

# save the output
cmd.png(
    "SecB_dG_render.png",
    width="10cm",
    dpi=300,
    ray=1,
)


# %%
# rotate for a different view
cmd.rotate("y", 90)
cmd.ipython_image()

# %%
