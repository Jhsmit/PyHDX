from pyhdx.datasets import HDXDataSet
from pyhdx.models import HDXMeasurement, HDXMeasurementSet
import numpy as np
from pathlib import Path
import yaml

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"

np.random.seed(43)


def test_load_from_yaml():
    yaml_pth = Path(input_dir / "data_states.yaml")
    hdx_spec = yaml.safe_load(yaml_pth.read_text())

    dataset = HDXDataSet.from_spec(hdx_spec, data_dir=input_dir)

    hdxm = HDXMeasurement.from_dataset(dataset, state="SecB_tetramer")
    assert isinstance(hdxm, HDXMeasurement)

    assert (
        hdxm.temperature
        == hdx_spec["states"]["SecB_tetramer"]["metadata"]["temperature"]["value"] + 273.15
    )

    assert hdxm.name == "SecB_tetramer"
    assert hdxm.state == "SecB WT apo"

    hdxm_set = HDXMeasurementSet.from_dataset(dataset)
    assert isinstance(hdxm_set, HDXMeasurementSet)
    assert hdxm_set.names == list(hdx_spec["states"].keys())
