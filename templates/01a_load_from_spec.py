"""Load a HDX-MS dataset from a .yaml spec file"""
# %%

from pathlib import Path

import numpy as np
import yaml

from pyhdx.batch_processing import StateParser

# %%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
np.random.seed(43)

yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
state_spec = yaml.safe_load(yaml_stream)

input_dir = current_dir.parent / "tests" / "test_data" / "input"

#%%

parser = StateParser(state_spec, input_dir)

#%%

hdxm = parser.load_hdxm('SecB_tetramer')
print(hdxm.timepoints)

#%%

hdxm = parser.load_hdxm('SecB_dimer')
print(hdxm)