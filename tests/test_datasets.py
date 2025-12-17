from pathlib import Path

import numpy as np

from pyhdx.models import HDXMeasurement, HDXMeasurementSet
from hdxms_datasets import load_dataset

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"

np.random.seed(43)

DATASET_ID = "HDX_D9096080"


def test_load_from_openHDX():
    dataset = load_dataset(input_dir / DATASET_ID)

    state = dataset.get_state("Tetramer")
    hdxm = HDXMeasurement.from_dataset(state)

    assert isinstance(hdxm, HDXMeasurement)

    assert hdxm.temperature == state.peptides[0].temperature

    assert hdxm.name == "Tetramer"
    assert hdxm.state == "SecB WT apo"

    hdxm_set = HDXMeasurementSet.from_dataset(dataset.states)
    assert isinstance(hdxm_set, HDXMeasurementSet)
    assert hdxm_set.names == [state.name for state in dataset.states]
