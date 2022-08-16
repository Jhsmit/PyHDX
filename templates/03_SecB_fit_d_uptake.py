# %%

from __future__ import annotations

from pathlib import Path

import numpy as np
import yaml

from pyhdx.batch_processing import StateParser
from pyhdx.fitting import fit_d_uptake, DUptakeFitResultSet

#%%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
np.random.seed(43)

yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
state_spec = yaml.safe_load(yaml_stream)

input_dir = current_dir.parent / "tests" / "test_data" / "input"

#%%
filters = [lambda df: df.query('exposure < 40.')]
parser = StateParser(state_spec, input_dir, data_filters=filters)


#%%
# load all hdx measurements, fit individually the D-uptake and then combine into one fit result
hdxm_set = parser.load_hdxmset()

results = []
for hdxm in hdxm_set:
    result = fit_d_uptake(hdxm, repeats=2)
    results.append(result)

results

#%%

fr = DUptakeFitResultSet(results)
fr.output

