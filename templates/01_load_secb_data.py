"""Load a HDX-MS dataset (peptide list with D-uptake per peptide, csv format)"""

# %%
from pathlib import Path

import numpy as np

from pyhdx.legacy import filter_peptides
from pyhdx.fileIO import read_dynamx
from pyhdx.models import HDXMeasurement
from pyhdx.process import apply_control, correct_d_uptake

# %%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
np.random.seed(43)

# %%

fpath = current_dir.parent / "tests" / "test_data" / "input" / "ecSecB_apo.csv"

# Load the full .csv file
df = read_dynamx(fpath)

fd = {"state": "Full deuteration control", "exposure": {"value": 0.167, "unit": "min"}}

# Filter out peptides for the full deuteration control
fd_df = filter_peptides(df, **fd)

# filter out peptides for the experiment with state "SecB WT apo"
peptides = filter_peptides(df, state="SecB WT apo")

# Apply FD control, returns only peptides in both FD control and experiment
peptides_control = apply_control(peptides, fd_df)

# Correct for N-terminal back exchanging residues and deuterium percentage of the buffer
peptides_corrected = correct_d_uptake(peptides_control, drop_first=1, d_percentage=90.0)

sequence = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"
temperature, pH = 273.15 + 30, 8.0

# %%

# Create HDX Measurement object with addtional experimental metadata (sequence, pH, temperature)
hdxm = HDXMeasurement(peptides_corrected, sequence=sequence, pH=pH, temperature=temperature)

# %%
hdxm.data
# %%
peptides_corrected
# %%
