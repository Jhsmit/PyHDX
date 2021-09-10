from pyhdx.batch_processing import yaml_to_hdxm, yaml_to_hdxmset
from pyhdx.models import HDXMeasurement, HDXMeasurementSet
import numpy as np
from pathlib import Path
import yaml

cwd = Path(__file__).parent
input_dir = cwd / 'test_data' / 'input'
output_dir = cwd / 'test_data' / 'output'

np.random.seed(43)


class TestBatchProcessing(object):

    def test_load_from_yaml(self):
        yaml_pth = Path(input_dir / 'data_states.yaml')
        data_dict = yaml.safe_load(yaml_pth.read_text())

        hdxm = yaml_to_hdxm(data_dict['SecB_tetramer'], data_dir=input_dir)
        assert isinstance(hdxm, HDXMeasurement)

        assert hdxm.metadata['temperature'] == data_dict['SecB_tetramer']['temperature']['value'] + 273.15
        assert hdxm.name == 'SecB WT apo'

        hdxm_set = yaml_to_hdxmset(data_dict, data_dir=input_dir)
        assert isinstance(hdxm_set, HDXMeasurementSet)
        assert hdxm_set.names == list(data_dict.keys())


