"""Load a HDX-MS dataset from an openHDX dataset"""

# %%
# %%
from pathlib import Path

from hdxms_datasets import load_dataset

from pyhdx.models import HDXMeasurement, HDXMeasurementSet

# %%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)

dataset_dir = current_dir.parent / "tests" / "test_data" / "input" / "HDX_D9096080"
dataset = load_dataset(dataset_dir)

dataset.states

# %%
# Load an HDX measurement by state name
hdxm = HDXMeasurement.from_dataset(dataset.get_state("Tetramer"))
print(hdxm)
print(hdxm.timepoints)


# %%
# Load an HDX measurement set from all states in the dataset
hdxm_set = HDXMeasurementSet.from_dataset(dataset.states)
