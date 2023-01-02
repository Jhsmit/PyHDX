"""Load a HDX-MS dataset from a .yaml HDX specification file"""
# %%
from pathlib import Path
import yaml

from pyhdx.batch_processing import StateParser

# %%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)

yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
hdx_spec = yaml.safe_load(yaml_stream)

input_dir = current_dir.parent / "tests" / "test_data" / "input"

#%%

parser = StateParser(hdx_spec, input_dir)

#%%
# Load an HDX measurement by state name
hdxm = parser.load_hdxm("SecB_tetramer")
print(hdxm.timepoints)

#%%

# Load an HDX measurement by index of states in the spec
hdxm = parser.load_hdxm(1)
print(hdxm)

#%%

# Load all HDX measurements as a HDX measurment set.
hdxm_set = parser.load_hdxmset()

print(hdxm_set)
