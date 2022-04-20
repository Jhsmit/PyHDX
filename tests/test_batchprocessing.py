from pyhdx.batch_processing import StateParser, JobParser
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
        data_dict = yaml.safe_load(yaml_pth.read_text())

        parser = StateParser(data_dict, data_src=input_dir)

        hdxm = parser.load_hdxm("SecB_tetramer")
        assert isinstance(hdxm, HDXMeasurement)

        assert (
            hdxm.metadata["temperature"]
            == data_dict["SecB_tetramer"]["temperature"]["value"] + 273.15
        )
        assert hdxm.name == "SecB WT apo"

        hdxm_set = parser.load_hdxmset()
        assert isinstance(hdxm_set, HDXMeasurementSet)
        assert hdxm_set.names == list(data_dict.keys())

    def test_load_job_parser(self):
        fit_output_dir = input_dir / "fit_result_output_1"
        if fit_output_dir.exists():
            shutil.rmtree(fit_output_dir, ignore_errors=True)

        job_spec = yaml.safe_load((input_dir / "jobfile.yaml").read_text())
        parser = JobParser(job_spec, cwd=input_dir)
        parser.execute()

        assert fit_output_dir.exists()
        shutil.rmtree(fit_output_dir, ignore_errors=True)
