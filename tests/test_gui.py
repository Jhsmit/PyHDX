from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries
from pyhdx.fileIO import txt_to_protein, txt_to_np
from pyhdx.panel.apps import _main_app, _diff_app
from pyhdx.panel.data_sources import DataSource
from pathlib import Path

import numpy as np

directory = Path(__file__).parent


TEST_PML = """set_color color_#0a0ac2, [10,10,194]
set_color color_#8c8c8c, [140,140,140]
color color_#0a0ac2, resi 10-17 + resi 19-25 + resi 27-28 + resi 30-37 + resi 39-57 + resi 62-84 + resi 86-94 + resi 100-102 + resi 104-107 + resi 109-113 + resi 115-123 + resi 125-129 + resi 131-133 + resi 138-155
color color_#8c8c8c, resi 18-18 + resi 26-26 + resi 29-29 + resi 38-38 + resi 58-61 + resi 85-85 + resi 95-99 + resi 103-103 + resi 108-108 + resi 114-114 + resi 124-124 + resi 130-130 + resi 134-137
"""


class TestMainGUISecB(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = directory / 'test_data' / 'ecSecB_apo.csv'
        cls.pmt = PeptideMasterTable(read_dynamx(cls.fpath))

        cls.state = 'SecB WT apo'
        cls.control = ('Full deuteration control', 0.167)
        cls.pmt.set_control(cls.control)

        states = cls.pmt.groupby_state()
        cls.series = states[cls.state]

        cls.prot_fit_result = txt_to_protein(directory / 'test_data' / 'ecSecB_torch_fit.txt')

        cls.ds_fit = DataSource(cls.prot_fit_result, name='global_fit', x='r_number', tags=['mapping', 'pfact', 'deltaG'],
                                renderer='circle', size=10)

    def test_load_files(self):
        with open(self.fpath, 'rb') as f:
            binary = f.read()

        tmpl, ctrl = _main_app()
        file_input_control = ctrl.control_panels['PeptideFileInputControl']

        file_input_control.file_selectors[0].value = binary
        file_input_control._action_load()
        assert file_input_control.fd_state == 'Full deuteration control'
        assert file_input_control.fd_exposure == 0.0

        file_input_control.fd_state = self.control[0]
        file_input_control.fdexposure = self.control[1]

        file_input_control.exp_state = self.state
        assert file_input_control.exp_exposures == [0.0, 0.167, 0.5, 1.0, 5.0, 10.0, 100.000008]
        file_input_control._action_parse()

        assert isinstance(ctrl.series, KineticsSeries)
        assert isinstance(ctrl.peptides, PeptideMasterTable)

    def test_coverage(self):
        tmpl, ctrl = _main_app()
        ctrl.series = self.series

        cov_figure = ctrl.figure_panels['CoverageFigure']
        renderer = cov_figure.figure.renderers[0]

        assert renderer.data_source.name == f'coverage_{self.series.state}'

    def test_initial_guesses_and_fit(self):
        tmpl, ctrl = _main_app()
        ctrl.series = self.series

        ctrl.control_panels['InitialGuessControl']._action_fit()
        assert 'half-life' in ctrl.sources.keys()

        fit_control = ctrl.control_panels['FitControl']
        fit_control.epochs = 10
        fit_control._do_fitting()

        deltaG_figure = ctrl.figure_panels['DeltaGFigure']
        renderer = deltaG_figure.figure.renderers[0]
        assert renderer.data_source.name == 'global_fit'

    def test_file_download_output(self):
        tmpl, ctrl = _main_app()
        ctrl.series = self.series
        ctrl.publish_data('global_fit', self.ds_fit)

        f_export = ctrl.control_panels['FileExportControl']
        f_export.target = 'global_fit'
        pml = f_export._make_pml('sadf')

        assert pml == TEST_PML

        io = f_export.linear_export_callback()
        val = io.getvalue()

        assert val.count('\n') == 148
        array = txt_to_np(io)
        assert np.allclose(f_export.export_dict['deltaG'], array['deltaG'], equal_nan=True)


class TestMainGUISimulated(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = directory / 'test_data' / 'simulated_data_uptake.csv'
        cls.pmt = PeptideMasterTable(read_dynamx(cls.fpath))

        cls.state = 'state1'
        cls.pmt.set_backexchange(0.)

        states = cls.pmt.groupby_state()
        cls.series = states[cls.state]

    def test_load_files(self):
        with open(self.fpath, 'rb') as f:
            binary = f.read()

        tmpl, ctrl = _main_app()
        file_input_control = ctrl.control_panels['PeptideFileInputControl']

        file_input_control.file_selectors[0].value = binary
        file_input_control._action_load()
        file_input_control.norm_mode = 'Theory'
        file_input_control.be_percent = 0.

        file_input_control.exp_state = self.state
        assert file_input_control.exp_exposures == [0.167, 0.5, 1.0, 5.0, 10.0, 30., 100.]
        file_input_control._action_parse()

        assert isinstance(ctrl.series, KineticsSeries)
        assert isinstance(ctrl.peptides, PeptideMasterTable)

    def test_coverage(self):
        tmpl, ctrl = _main_app()
        ctrl.series = self.series

    def test_initial_guesses(self):
        tmpl, ctrl = _main_app()
        ctrl.cluster = None
        ctrl.series = self.series

        initial_guess = ctrl.control_panels['InitialGuessControl']
        initial_guess._action_fit()
        assert 'half-life' in ctrl.sources.keys()

        initial_guess.fitting_model = 'Association'
        initial_guess._action_fit()
        assert 'fit1' in ctrl.sources.keys()

    def test_global_fit(self):
        tmpl, ctrl = _main_app()
        ctrl.cluster = None
        ctrl.series = self.series

        initial_guess = ctrl.control_panels['InitialGuessControl']
        initial_guess._action_fit()

        fit = ctrl.control_panels['FitControl']
        fit.epochs = 20
        fit._do_fitting()
        assert 'global_fit' in ctrl.sources.keys()


class TestDiffApp(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = directory / 'test_data' / 'ecSecB_torch_fit.txt'

        with open(cls.fpath, 'rb') as f_obj:
            cls.file_binary = f_obj.read()

    def test_app(self):
        tmpl, ctrl = _diff_app()

        f_input = ctrl.control_panels['MappingFileInputControl']
        f_input._widget_dict['input_file'].filename = str(self.fpath)
        f_input.input_file = self.file_binary
        assert f_input.dataset_name == 'ecSecB_torch_fit'
        f_input.dataset_name = 'DS1'
        f_input._action_add_dataset()
        assert f_input.input_file == b''
        assert f_input.dataset_name == ''

        f_input = ctrl.control_panels['MappingFileInputControl']
        f_input._widget_dict['input_file'].filename = str(self.fpath)
        f_input.input_file = self.file_binary
        f_input.dataset_name = 'DS2'
        f_input._action_add_dataset()

        diff = ctrl.control_panels['DifferenceControl']
        diff.dataset_1 = 'DS1'
        diff.dataset_2 = 'DS2'

        comparison_name = 'Diff_ds1_ds2'
        diff.comparison_name = comparison_name
        quantity_objects = diff.param['comparison_quantity'].objects
        assert quantity_objects == sorted(['_deltaG', 'deltaG', 'pfact'])

        diff.comparison_quantity = 'deltaG'
        diff._action_add_comparison()
        assert diff.comparison_name == ''
        assert comparison_name in ctrl.sources.keys()
        assert diff.comparison_list is None

