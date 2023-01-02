from pyhdx.batch_processing import StateParser
from pyhdx.models import HDXMeasurement, HDXMeasurementSet
import numpy as np
from pathlib import Path
import yaml
import shutil

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"

np.random.seed(43)


class TestBatchProcessing(object):
    def test_load_from_yaml(self):
        yaml_pth = Path(input_dir / "data_states.yaml")
        hdx_spec = yaml.safe_load(yaml_pth.read_text())

        parser = StateParser(hdx_spec, data_src=input_dir)

        hdxm = parser.load_hdxm("SecB_tetramer")
        assert isinstance(hdxm, HDXMeasurement)

        assert (
            hdxm.temperature
            == hdx_spec["states"]["SecB_tetramer"]["metadata"]["temperature"]["value"]
            + 273.15
        )

        assert hdxm.name == "SecB_tetramer"
        assert hdxm.state == "SecB WT apo"

        hdxm_set = parser.load_hdxmset()
        assert isinstance(hdxm_set, HDXMeasurementSet)
        assert hdxm_set.names == list(hdx_spec["states"].keys())
