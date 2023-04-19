from pathlib import Path

from pyhdx.fileIO import csv_to_dataframe
from pyhdx.config import cfg
from pyhdx.plot import (
    dG_scatter_figure,
    linear_bars_figure,
    rainbowclouds_figure,
    ddG_scatter_figure,
)
import proplot as pplt

# %%
cwd = Path(__file__).parent
root_dir = cwd.parent
web_data_dir = root_dir / "tests" / "test_data" / "output" / "web"


# %%
dG_df = csv_to_dataframe(web_data_dir / "dG.csv")
dG_df
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
pplt.show()

# %%

linear_bars_figure(plot_data)
pplt.show()

# %%

rainbowclouds_figure(plot_data)
pplt.show()

# %%
protein_states = plot_data.columns.get_level_values(0).unique()
ddG_scatter_figure(plot_data, reference=protein_states[0])
