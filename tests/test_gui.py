from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries
from pyhdx.panel.apps import _main_app, _diff_app
from pathlib import Path

directory = Path(__file__).parent


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

    def test_load_files(self):
        with open(self.fpath, 'rb') as f:
            binary = f.read()

        tmpl, ctrl = _main_app()
        file_input_control = ctrl.control_panels['PeptideFileInputControl']

        file_input_control.file_selectors[0].value = binary
        file_input_control._action_load()
        assert file_input_control.norm_state == 'Full deuteration control'
        assert file_input_control.norm_exposure == 0.0

        file_input_control.norm_state = self.control[0]
        file_input_control.norm_exposure = self.control[1]

        file_input_control.exp_state = self.state
        assert file_input_control.exp_exposures == [0.0, 0.167, 0.5, 1.0, 5.0, 10.0, 100.000008]
        file_input_control._action_parse()

        assert isinstance(ctrl.series, KineticsSeries)
        assert isinstance(ctrl.peptides, PeptideMasterTable)

    def test_coverage(self):
        tmpl, ctrl = _main_app()
        ctrl.series = self.series

    def test_initial_guesses(self):
        tmpl, ctrl = _main_app()
        ctrl.series = self.series

        ctrl.control_panels['InitialGuessControl']._action_fit()
        assert 'half-life' in ctrl.sources.keys()


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

    def test_tf_fit(self):
        tmpl, ctrl = _main_app()
        ctrl.cluster = None
        ctrl.series = self.series

        initial_guess = ctrl.control_panels['InitialGuessControl']
        initial_guess._action_fit()

        fit = ctrl.control_panels['FitControl']
        fit._do_fitting()
        assert 'pfact' in ctrl.sources.keys()


class TestDiffApp(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = directory / 'test_data' / 'SecB WT apo_pfact_linear.txt'

        with open(cls.fpath, 'rb') as f_obj:
            cls.file_binary = f_obj.read()

    def test_app(self):
        tmpl, ctrl = _diff_app()

        f_input = ctrl.control_panels['MappingFileInputControl']
        f_input._widget_dict['input_file'].filename = str(self.fpath)
        f_input.input_file = self.file_binary
        assert f_input.dataset_name == 'SecB WT apo_pfact_linear'
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
        assert quantity_objects == sorted(['log_P_full', 'log_P', 'deltaG'])

        diff.comparison_quantity = 'deltaG'
        diff._action_add_comparison()
        assert diff.comparison_name == ''
        assert comparison_name in ctrl.sources.keys()
        assert diff.comparison_list is None

