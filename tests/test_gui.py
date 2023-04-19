from pathlib import Path

import numpy as np
import pytest
import torch
import yaml
from distributed.utils_test import cluster
import pandas as pd

from pyhdx import read_dynamx, HDXMeasurement
from pyhdx.config import cfg
from pyhdx.fileIO import csv_to_dataframe
from pyhdx.process import filter_peptides, apply_control, correct_d_uptake
from pyhdx.support import hash_dataframe, hash_array
from pyhdx.web.apps import main_app, rfu_app
from pyhdx.web.utils import load_state_rfu

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"


TEST_PML = """set_color color_#0a0ac2, [10,10,194]
set_color color_#8c8c8c, [140,140,140]
color color_#0a0ac2, resi 10-17 + resi 19-25 + resi 27-28 + resi 30-37 + resi 39-57 + resi 62-84 + resi 86-94 + resi 100-102 + resi 104-107 + resi 109-113 + resi 115-123 + resi 125-129 + resi 131-133 + resi 138-155
color color_#8c8c8c, resi 1-9 + resi 18-18 + resi 26-26 + resi 29-29 + resi 38-38 + resi 58-61 + resi 85-85 + resi 95-99 + resi 103-103 + resi 108-108 + resi 114-114 + resi 124-124 + resi 130-130 + resi 134-137
"""

test_port = 55432
np.random.seed(43)
torch.manual_seed(43)

@pytest.fixture
def ppix_spec() -> dict:
    return yaml.safe_load(Path(input_dir / "PpiX_states.yaml").read_text())

@pytest.fixture
def secb_spec() -> dict:
    return yaml.safe_load(Path(input_dir / "data_states.yaml").read_text())


class TestMainGUISecB(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = input_dir / "ecSecB_apo.csv"
        df = read_dynamx(cls.fpath)

        fd = {
            "state": "Full deuteration control",
            "exposure": {"value": 0.167, "unit": "min"},
        }

        fd_df = filter_peptides(df, **fd)
        peptides = filter_peptides(df, state="SecB WT apo")  # , query=["exposure != 0."])
        peptides_control = apply_control(peptides, fd_df)
        peptides_corrected = correct_d_uptake(peptides_control)

        cls.temperature, cls.pH = 273.15 + 30, 8.0
        cls.hdxm = HDXMeasurement(
            peptides_corrected, temperature=cls.temperature, pH=cls.pH, c_term=155
        )

    def test_load_single_file(self):
        with open(self.fpath, "rb") as f:
            binary = f.read()

        ctrl, tmpl = main_app()
        src = ctrl.sources["main"]
        input_control = ctrl.control_panels["PeptideRFUFileInputControl"]

        input_control.widgets["input_files"].filename = ["ecSecB_apo.csv"]
        input_control.input_files = [binary]
        assert input_control.fd_state == "Full deuteration control"
        assert input_control.fd_exposure == 0.0

        input_control.fd_state = "Full deuteration control"
        input_control.fd_exposure = 10.020000000000001

        input_control.exp_state = "SecB WT apo"
        timepoints = list(np.array([0.167, 0.5, 1.0, 5.0, 10.0, 100.000008]) * 60)
        assert input_control.exp_exposures == timepoints
        input_control._add_single_dataset_spec()
        input_control._action_load_datasets()

        assert "SecB WT apo" in src.hdxm_objects
        hdxm = src.hdxm_objects["SecB WT apo"]

        assert hdxm.Nt == 6
        assert hdxm.Np == 63
        assert hdxm.Nr == 145

        assert np.nanmean(hdxm.rfu_residues) == pytest.approx(0.6335831166442542)

    def test_batch_input(self):
        filenames = ["ecSecB_apo.csv", "ecSecB_dimer.csv"]
        file_dict = {fname: (input_dir / fname).read_bytes() for fname in filenames}

        ctrl, tmpl = main_app()

        input_control = ctrl.control_panels["PeptideRFUFileInputControl"]
        input_control.input_mode = "Batch"
        input_control.widgets["input_files"].filename = list(file_dict.keys())
        input_control.input_files = list(file_dict.values())

        input_control.batch_file = Path(input_dir / "data_states.yaml").read_bytes()

        input_control._action_load_datasets()

        src = ctrl.sources["main"]
        assert len(src.hdxm_objects) == 2
        # ... additional tests

    # @pytest.mark.skip(reason="Fails in GitHub Actions")
    def test_web_fitting(self):
        filenames = ["ecSecB_apo.csv", "ecSecB_dimer.csv"]
        file_dict = {fname: (input_dir / fname).read_bytes() for fname in filenames}

        ctrl, tmpl = main_app()

        input_control = ctrl.control_panels["PeptideRFUFileInputControl"]
        # input_control.input_mode = "Batch"
        input_control.widgets["input_files"].filename = list(file_dict.keys())
        input_control.input_files = list(file_dict.values())

        filenames = ["ecSecB_apo.csv", "ecSecB_dimer.csv"]
        input_control.widgets["input_files"].filename = filenames

        input_control.fd_state = "Full deuteration control"
        input_control.fd_exposure = 0.167 * 60

        input_control.exp_state = "SecB WT apo"
        input_control.measurement_name = "testname_123"
        input_control._add_single_dataset_spec()

        input_control.exp_file = "ecSecB_dimer.csv"
        input_control.exp_state = "SecB his dimer apo"
        input_control.measurement_name = "SecB his dimer apo"  # todo catch error duplicate name
        input_control._add_single_dataset_spec()

        input_control._action_load_datasets()

        assert "testname_123" in ctrl.sources["main"].hdxm_objects.keys()
        assert "SecB his dimer apo" in ctrl.sources["main"].hdxm_objects.keys()

        rfu_df = ctrl.sources["main"].get_table("rfu")
        assert rfu_df.shape == (145, 24)
        assert rfu_df.columns.nlevels == 3

        initial_guess = ctrl.control_panels["InitialGuessControl"]
        initial_guess._action_fit()


def test_web_load(secb_spec):
    ctrl, tmpl = main_app()

    file_input = ctrl.control_panels["PeptideRFUFileInputControl"]
    states = ['SecB_tetramer', 'SecB_dimer']
    load_state_rfu(file_input, secb_spec, data_dir=input_dir, states=states)

    file_input._action_load_datasets()
    assert len(file_input.src.hdxm_objects) == 2

    file_export = ctrl.control_panels["FileExportControl"]

    # check table output
    file_export.table = 'peptides'
    sio = file_export.table_export_callback()
    sio.seek(0)
    lines_test = sio.read().split("\n")
    lines_ref = (output_dir / 'main_web' / 'peptides.csv').read_text().split("\n")

    for lt, lr in zip(lines_test[2:], lines_ref[2:]):
        assert lt == lr

    # check rfu table output
    file_export.table = 'rfu'
    sio = file_export.table_export_callback()
    sio.seek(0)
    lines_test = sio.read().split("\n")
    lines_ref = (output_dir / 'main_web' / 'rfu.csv').read_text().split("\n")

    for lt, lr in zip(lines_test[2:], lines_ref[2:]):
        assert lt == lr

    #check color table output
    sio = file_export.color_export_callback()
    sio.seek(0)
    lines_test = sio.read().split("\n")
    lines_ref = (output_dir / 'main_web' / 'rfu_colors.csv').read_text().split("\n")

    for lt, lr in zip(lines_test[2:], lines_ref[2:]):
        assert lt == lr

def test_rfu(ppix_spec):
    """Test the RFU app"""
    ctrl, tmpl = rfu_app()

    file_input = ctrl.control_panels["PeptideRFUFileInputControl"]
    states = ['PpiA_Folding', 'PpiB_Folding']
    load_state_rfu(file_input, ppix_spec, data_dir = input_dir, states=states)

    file_input._action_load_datasets()
    assert len(file_input.src.hdxm_objects) == 2

    file_export = ctrl.control_panels["FileExportControl"]

    # check table output
    file_export.table = 'peptides'
    sio = file_export.table_export_callback()
    sio.seek(0)
    lines_test = sio.read().split("\n")
    lines_ref = (output_dir / 'rfu_web' / 'peptides.csv').read_text().split("\n")

    for lt, lr in zip(lines_test[2:], lines_ref[2:]):
        assert lt == lr

    # check rfu table output
    file_export.table = 'rfu'
    sio = file_export.table_export_callback()
    sio.seek(0)
    lines_test = sio.read().split("\n")
    lines_ref = (output_dir / 'rfu_web' / 'rfu.csv').read_text().split("\n")

    for lt, lr in zip(lines_test[2:], lines_ref[2:]):
        assert lt == lr




        # with cluster() as (s, [a, b]):
        #     conf.set("cluster", "scheduler_address", s["address"])
        #
        #     fit_control = ctrl.control_panels["FitControl"]
        #     fit_control.fit_mode = "Batch"
        #     fit_control.epochs = 10
        #
        #     fit_control.fit_name = "testfit_1"
        #     fit_control._action_fit()
        #
        #     table_names = [
        #         "loss",
        #         "peptide_mse",
        #         "rates",
        #         "dG_fits",
        #         "rfu_residues",
        #         "peptides",
        #         "d_calc",
        #     ]
        #     for name in table_names:
        #         ref = csv_to_dataframe(output_dir / "gui" / f"{name}.csv")
        #         test = ctrl.sources["main"].get_table(name)
        #
        #         for test_col, ref_col in zip(test.columns, ref.columns):
        #             test_values = test[test_col].to_numpy()
        #             if np.issubdtype(test_values.dtype, np.number):
        #                 assert np.allclose(
        #                     test_values, ref[ref_col].to_numpy(), equal_nan=True
        #                 )
        #
        #     fit_control.guess_mode = "One-to-many"
        #     fit_control.fit_name = "testfit_2"
        #     fit_control._action_fit()
        #
        #     fit_control.fit_mode = "Single"
        #     fit_control.guess_mode = "One-to-one"
        #     fit_control.fit_name = "testfit_3"
        #     fit_control._action_fit()
        #
        #     fit_control.initial_guess = "testfit_2"
        #     fit_control.guess_mode = "One-to-many"
        #     fit_control.guess_state = "testname_123"
        #     fit_control.fit_name = "testfit_4"
        #     fit_control._action_fit()
        #
        # color_transform_control = ctrl.control_panels["ColorTransformControl"]
        # color_transform_control._action_otsu()
        #
        # fit_result = ctrl.sources["main"].get_table("dG_fits")
        # values = fit_result["testfit_1"]["testname_123"]["dG"]
        # color_transform_control.quantity = "dG"
        # cmap, norm = color_transform_control.get_cmap_and_norm()
        # colors = cmap(norm(values), bytes=True)
        #
        # h = hash_array(colors, method="md5")
        # assert h == "bac3602f877a53abd94be8bb5f9b72ec"
        #
        # color_transform_control.mode = "Continuous"
        # color_transform_control._action_linear()
        #
        # cmap, norm = color_transform_control.get_cmap_and_norm()
        # colors = cmap(norm(values), bytes=True)
        #
        # h = hash_array(colors, method="md5")
        #
        # assert h == "123085ba16b3a9374595b734f9e675e6"
        #
        # value_widget = color_transform_control.widgets["value_1"]
        # value_widget.value = 25
        # cmap, norm = color_transform_control.get_cmap_and_norm()
        # assert norm.vmin == 25
        #
        # color_transform_control.mode = "Colormap"
        # color_transform_control.library = "colorcet"
        # color_transform_control.colormap = "CET_C1"
        # cmap, norm = color_transform_control.get_cmap_and_norm()
        #
        # colors = cmap(norm(values), bytes=True)
        #
        # h = hash_array(colors, method="md5")
        # assert h == "0628b46e7975ed57490e84c169bc81ad"

        # Future tests: check if renderers are present in figures
        # cov_figure = ctrl.figure_panels['CoverageFigure']
        # renderer = cov_figure.figure.renderers[0]
        # assert renderer.data_source.name == f'coverage_{self.series.state}'