"""Load HDX-MS data from yaml spec and perform initial guess of exchange rates"""
from pyhdx import HDXMeasurement
from pyhdx.datasets import HDXDataSet
from pathlib import Path
from pyhdx.fitting import fit_rates_weighted_average
import yaml
from pyhdx.local_cluster import default_client
from pyhdx.fileIO import dataframe_to_file

current_dir = Path(__file__).parent
output_dir = current_dir / "guesses"
output_dir.mkdir(exist_ok=True)
data_dir = current_dir.parent / "tests" / "test_data" / "input"
yaml_stream = Path(current_dir / "yaml_files" / "SecB.yaml").read_text()
hdx_spec = yaml.safe_load(yaml_stream)

dataset = HDXDataSet.from_spec(
    hdx_spec,
    data_dir=data_dir,
)

# Requires local_cluster.py to be running (or other Dask client on default address in config)
client = default_client()

for name in dataset.states:
    print(name)
    hdxm = HDXMeasurement.from_dataset(dataset, state=name)

    # Save sequence info + intrinsic rates
    dataframe_to_file(
        output_dir / f"{name}_intrinsic_rates.csv", hdxm.coverage.protein, fmt="pprint"
    )

    fr = fit_rates_weighted_average(hdxm, client=client)
    dataframe_to_file(output_dir / f"{name}_rates_guess.csv", fr.output)
