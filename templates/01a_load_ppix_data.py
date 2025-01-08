"""Load a HDX-MS dataset with both a FD and ND control"""

from pathlib import Path

import numpy as np

from pyhdx import read_dynamx
from pyhdx.process import apply_control, filter_peptides

# %%

root_dir = Path(__file__).parent.parent
np.random.seed(43)

fpath = root_dir / "tests" / "test_data" / "input" / "PpiA_folding.csv"

# Load the full .csv file
df = read_dynamx(fpath)
df
# %%

fd_spec = {"state": "Native", "exposure": {"value": 86400, "unit": "s"}}
nd_spec = {"state": "FD", "exposure": {"value": 0.6, "unit": "s"}}

# Filter out peptides for the full deuteration control
fd_df = filter_peptides(df, **fd_spec)

# Filter out peptides for the non-deuterated control
nd_df = filter_peptides(df, **nd_spec)

# Filter out peptides for the experiment with state "Folding"
peptides = filter_peptides(df, state="Folding")

# %%
# Apply the FD and ND controls. This adds the columns 'rfu' and 'rfu_sd' (Relative Fractional Uptake)
# as well as 'fd_uptake' and 'nd_uptake', and their standard deviations.
peptides_control = apply_control(peptides, fd_control=fd_df, nd_control=nd_df)
peptides_control.columns

# %%
