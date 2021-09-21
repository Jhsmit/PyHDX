from pyhdx import PeptideMasterTable, read_dynamx, HDXMeasurement
from pyhdx.fileIO import csv_to_protein
from pyhdx.web.apps import main_app#, diff_app
from pyhdx.config import ConfigurationSettings
from pathlib import Path
import torch
import numpy as np
import pytest
import time

cwd = Path(__file__).parent
input_dir = cwd / 'test_data' / 'input'
output_dir = cwd / 'test_data' / 'output'


TEST_PML = """set_color color_#0a0ac2, [10,10,194]
set_color color_#8c8c8c, [140,140,140]
color color_#0a0ac2, resi 10-17 + resi 19-25 + resi 27-28 + resi 30-37 + resi 39-57 + resi 62-84 + resi 86-94 + resi 100-102 + resi 104-107 + resi 109-113 + resi 115-123 + resi 125-129 + resi 131-133 + resi 138-155
color color_#8c8c8c, resi 1-9 + resi 18-18 + resi 26-26 + resi 29-29 + resi 38-38 + resi 58-61 + resi 85-85 + resi 95-99 + resi 103-103 + resi 108-108 + resi 114-114 + resi 124-124 + resi 130-130 + resi 134-137
"""

test_port = 55432
np.random.seed(43)
torch.manual_seed(43)


class TestMainGUISecB(object):
    @classmethod
    def setup_class(cls):
        cls.fpath =input_dir / 'ecSecB_apo.csv'
        cls.pmt = PeptideMasterTable(read_dynamx(cls.fpath))

        cls.state = 'SecB WT apo'
        cls.control = ('Full deuteration control', 0.167*60)
        cls.pmt.set_control(cls.control)

        state_data = cls.pmt.get_state(cls.state)
        cls.temperature, cls.pH = 273.15 + 30, 8.
        cls.series = HDXMeasurement(state_data, temperature=cls.temperature, pH=cls.pH)
        cls.prot_fit_result = csv_to_protein(output_dir / 'ecSecB_torch_fit.csv')

        cfg = ConfigurationSettings()
        cfg.set('cluster', 'scheduler_address', f'127.0.0.1:{test_port}')
        #cfg.set('cluster', 'port', str(test_port))

    def test_load_single_file(self):
        with open(self.fpath, 'rb') as f:
            binary = f.read()

        ctrl = main_app()
        file_input_control = ctrl.control_panels['PeptideFileInputControl']

        file_input_control.input_files = [binary]
        assert file_input_control.fd_state == 'Full deuteration control'
        assert file_input_control.fd_exposure == 0.0

        file_input_control.fd_state = self.control[0]
        file_input_control.fd_exposure = self.control[1]

        file_input_control.exp_state = self.state
        timepoints = list(np.array([0.0, 0.167, 0.5, 1.0, 5.0, 10.0, 100.000008])*60)
        assert file_input_control.exp_exposures == timepoints
        file_input_control._action_add_dataset()

        assert self.state in ctrl.data_objects
        hdxm = ctrl.data_objects[self.state]

        assert hdxm.Nt == 7
        assert hdxm.Np == 63
        assert hdxm.Nr == 146

        file_input_control.input_files = [binary]
        assert file_input_control.fd_state == 'Full deuteration control'
        assert file_input_control.fd_state == self.control[0]
        assert file_input_control.fd_exposure == self.control[1]

        file_input_control.exp_state = self.state
        timepoints = list(np.array([0.0, 0.5, 1.0, 5.0, 100.000008])*60)
        file_input_control.exp_exposures = timepoints
        print('second add')
        #file_input_control.dataset_name = 'test123'
        # Add dataset with the name name
        file_input_control._action_add_dataset()
        assert len(ctrl.data_objects) == 1

        name = 'second_dataset'
        file_input_control.dataset_name = name
        file_input_control._action_add_dataset()

        assert name in ctrl.data_objects
        hdxm = ctrl.data_objects[name]

        assert hdxm.Nt == 5
        assert hdxm.Np == 63
        assert hdxm.Nr == 146

        assert np.nanmean(hdxm.rfu_residues) == pytest.approx(0.521232427211)

    @pytest.mark.skip(reason="Hangs in GitHub Actions")
    def test_batch_mode(self):
        fpath_1 = input_dir / 'ecSecB_apo.csv'
        fpath_2 = input_dir / 'ecSecB_dimer.csv'

        fpaths = [fpath_1, fpath_2]
        files = [p.read_bytes() for p in fpaths]

        ctrl = main_app(client=None)
        file_input = ctrl.control_panels['PeptideFileInputControl']

        file_input.input_files = files
        file_input.fd_state = 'Full deuteration control'
        file_input.fd_exposure = 0.167*60

        file_input.exp_state = 'SecB WT apo'
        file_input.dataset_name = 'testname_123'
        file_input._action_add_dataset()

        #assert ....

        file_input.exp_state = 'SecB his dimer apo'
        file_input.dataset_name = 'SecB his dimer apo'  # todo catch error duplicate name
        file_input._action_add_dataset()

        initial_guess = ctrl.control_panels['InitialGuessControl']
        initial_guess._action_fit()

        # Wait until fitting futures are completed
        while len(ctrl.future_queue) > 0:
            ctrl.check_futures()
            time.sleep(0.1)


        #assert ....
        guesses = ctrl.sources['dataframe'].get('rates')


        fit_control = ctrl.control_panels['FitControl']
        fit_control.epochs = 10

        fit_control.fit_name = 'testfit_1'
        fit_control._action_fit()

        # Wait until fitting futures are completed
        while len(ctrl.future_queue) > 0:
            ctrl.check_futures()
            time.sleep(0.1)
        #assert ....

        table = ctrl.sources['dataframe'].get('global_fit')

        # Test classification
        # todo test log space
        # todo probably values should be fixed otherwise tests are co-dependent
        values = table['testfit_1']['testname_123']['deltaG']
        classification = ctrl.control_panels['ClassificationControl']
        classification._action_otsu()

        cmap, norm = classification.get_cmap_and_norm()
        colors = cmap(norm(values), bytes=True)

        assert colors.sum() == 68474
        assert colors.std() == 109.39692637364178

        classification.mode = 'Continuous'
        classification._action_linear()

        cmap, norm = classification.get_cmap_and_norm()
        colors = cmap(norm(values), bytes=True)

        assert colors.sum() == 73090
        assert colors.std() == 89.41501408256475

        value_widget = classification.widgets['value_2']
        value_widget.value = 10e3

        cmap, norm = classification.get_cmap_and_norm()
        colors = cmap(norm(values), bytes=True)

        assert colors.sum() == 73097
        assert colors.std() == 91.9688274382922

        classification.mode = 'Color map'
        classification.library = 'colorcet'
        classification.color_map = 'CET_C1'
        cmap, norm = classification.get_cmap_and_norm()

        colors = cmap(norm(values), bytes=True)

        assert colors.sum() == 117289
        assert colors.std() == 64.90120978241222

        #
        # cov_figure = ctrl.figure_panels['CoverageFigure']
        # renderer = cov_figure.figure.renderers[0]
        #
        # assert renderer.data_source.name == f'coverage_{self.series.state}'

    def test_initial_guesses_and_fit(self):
        ctrl = main_app(client=None)
        ctrl.data_objects[self.series.name] = self.series

        ctrl.control_panels['InitialGuessControl']._action_fit()

        # todo Add tests
        # assert 'half-life' in ctrl.sources.keys()
        #
        # fit_control = ctrl.control_panels['FitControl']
        # fit_control.epochs = 10
        # fit_control._do_fitting()
        #
        # deltaG_figure = ctrl.figure_panels['DeltaGFigure']
        # renderer = deltaG_figure.figure.renderers[0]
        # assert renderer.data_source.name == 'global_fit'

    def test_file_download_output(self):
        ctrl = main_app()

        # todo Add tests
        # ctrl.series = self.series
        # ctrl.publish_data('global_fit', self.ds_fit)
        #
        # f_export = ctrl.control_panels['FileExportControl']
        # f_export.target = 'global_fit'
        # pml = f_export._make_pml('sadf')
        #
        # assert pml == TEST_PML
        #
        # io = f_export.linear_export_callback()
        # val = io.getvalue()
        #
        # assert val.count('\n') == 148


@pytest.mark.skip(reason="Diff app needs overhaul to lumen framework")
class TestDiffApp(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = directory / 'test_data' / 'ecSecB_torch_fit.csv'

        with open(cls.fpath, 'rb') as f_obj:
            cls.file_binary = f_obj.read()

    def test_app(self):
        ctrl = diff_app()

        f_input = ctrl.control_panels['MappingFileInputControl']
        f_input.widget_dict['input_file'].filename = str(self.fpath)
        f_input.input_file = self.file_binary
        assert f_input.dataset_name == 'ecSecB_torch_fit'
        f_input.dataset_name = 'DS1'
        f_input._action_add_dataset()
        assert f_input.input_file == b''
        assert f_input.dataset_name == ''

        f_input = ctrl.control_panels['MappingFileInputControl']
        f_input.widget_dict['input_file'].filename = str(self.fpath)
        f_input.input_file = self.file_binary
        f_input.dataset_name = 'DS2'
        f_input._action_add_dataset()

        diff = ctrl.control_panels['DifferenceControl']
        diff.dataset_1 = 'DS1'
        diff.dataset_2 = 'DS2'

        comparison_name = 'Diff_ds1_ds2'
        diff.comparison_name = comparison_name
        quantity_objects = diff.param['comparison_quantity'].objects
        assert quantity_objects == sorted(['_deltaG', 'covariance', 'deltaG', 'pfact'])

        diff.comparison_quantity = 'deltaG'
        diff._action_add_comparison()
        assert diff.comparison_name == ''
        assert comparison_name in ctrl.sources.keys()
        assert diff.comparison_list is None
#
