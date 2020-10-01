import os

from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries
from pyhdx.panel.apps import _main_app

from pyhdx.panel.log import get_default_handler
import sys
directory = os.path.dirname(__file__)


class TestGui(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = os.path.join(directory, 'test_data', 'ecSecB_apo.csv')
        cls.pmt = PeptideMasterTable(read_dynamx(cls.fpath))

        cls.state = 'SecB WT apo'
        cls.control = ('Full deuteration control', 0.167)
        cls.pmt.set_control(cls.control)

        states = cls.pmt.groupby_state()
        cls.series = states[cls.state]
        cls.series.make_uniform()

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