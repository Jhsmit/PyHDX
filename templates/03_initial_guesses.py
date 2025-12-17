"""Load HDX-MS data from an openHDX dataset and perform initial guess of exchange rates"""

from pyhdx import HDXMeasurementSet
from hdxms_datasets import load_dataset
from pathlib import Path
from pyhdx.fitting import fit_rates_weighted_average
from pyhdx.local_cluster import default_client
from pyhdx.fileIO import dataframe_to_file

current_dir = Path(__file__).parent
output_dir = current_dir / "guesses"
output_dir.mkdir(exist_ok=True)

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)

dataset_dir = current_dir.parent / "tests" / "test_data" / "input" / "HDX_D9096080"
dataset = load_dataset(dataset_dir)

# Requires local_cluster.py to be running (or other Dask client on default address in config)
client = default_client()

hdxm_set = HDXMeasurementSet.from_dataset(dataset.states)

for hdxm in hdxm_set:
    print(hdxm.name)

    # Save sequence info + intrinsic rates
    dataframe_to_file(
        output_dir / f"{hdxm.name}_intrinsic_rates.csv", hdxm.coverage.protein, fmt="pprint"
    )

    fr = fit_rates_weighted_average(hdxm, client=client)
    dataframe_to_file(output_dir / f"{hdxm.name}_rates_guess.csv", fr.output)
