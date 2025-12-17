"""Load a HDX-MS dataset from an openHDX dataset"""

# %%
# %%
from pathlib import Path

from hdxms_datasets import load_dataset

from pyhdx.fitting import fit_d_uptake
from pyhdx.models import HDXMeasurement, HDXMeasurementSet

# %%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)

dataset_dir = current_dir.parent / "tests" / "test_data" / "input" / "HDX_D9096080"
dataset = load_dataset(dataset_dir)

dataset.states


# %%

# zip_pth = Path(r"C:\Users\jhsmi\repos\mine\hdxms-datasets\tests\datasets\HDX_3BAE2080.zip")
# dataset = load_dataset(zip_pth)
# dataset

# %%
# Load an HDX measurement by state name
hdxm = HDXMeasurement.from_dataset(dataset.get_state(0))
hdxm = HDXMeasurement.from_dataset(dataset.get_state("Tetramer"), d_percentage=100.0, drop_first=1)
print(hdxm)
print(hdxm.timepoints)


# %%

fr = fit_d_uptake(hdxm, r1=0.5, repeats=3, verbose=False)
fr.output

# %%
from pyhdx.fileIO import dataframe_to_file

dataframe_to_file(output_dir / "hdxm_data.csv", fr.output)

# %%
# Load an HDX measurement set from all states in the dataset
hdxm_set = HDXMeasurementSet.from_dataset(dataset.states)
