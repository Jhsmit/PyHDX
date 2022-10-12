# %%

from pathlib import Path
import yaml

from pyhdx.batch_processing import StateParser
from pyhdx.fitting import fit_d_uptake, DUptakeFitResultSet

#%%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
hdx_spec = yaml.safe_load(yaml_stream)

#%%
# Modify the HDX spec to Limit the data to < 40. seconds exposure
for state, state_spec in hdx_spec['states'].items():
    peptide_spec = state_spec['peptides']
    peptide_spec['experiment']['query'] = ['exposure < 40.']

input_dir = current_dir.parent / "tests" / "test_data" / "input"
parser = StateParser(hdx_spec, input_dir)

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

