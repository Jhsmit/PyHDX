"""Load a HDX-MS dataset from a .yaml HDX specification file"""
# %%
from pathlib import Path
import yaml

from pyhdx.datasets import HDXDataSet

from pyhdx import HDXMeasurement, HDXMeasurementSet

# %%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)

yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
hdx_spec = yaml.safe_load(yaml_stream)
input_dir = current_dir.parent / "tests" / "test_data" / "input"

# %%

dataset = HDXDataSet.from_spec(
    hdx_spec,
    data_dir=input_dir,
)

print(dataset.describe())

# %%
# Load an HDX measurement by state name
hdxm = HDXMeasurement.from_dataset(dataset, state="SecB_tetramer")
print(hdxm)
print(hdxm.timepoints)

# # Load an HDX measurement by index of states in the spec

hdxm = HDXMeasurement.from_dataset(dataset, state=1)
print(str(hdxm))
print(hdxm.timepoints)


# %%
hdxm_set = HDXMeasurementSet.from_dataset(dataset)
print(hdxm_set)
