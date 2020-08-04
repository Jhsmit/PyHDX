import os

from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries
from pyhdx.panel.controller import Controller

directory = os.path.dirname(__file__)

class TestGui(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        cls.pmt = PeptideMasterTable(read_dynamx(cls.fpath))

        cls.state = 'PpiANative'
        cls.control = ('PpiA-FD', 0.167)
        cls.pmt.set_control(cls.control)

        states = cls.pmt.groupby_state()
        cls.series = states[cls.state]
        cls.series.make_uniform()

    def test_load_files(self):
        with open(self.fpath, 'rb') as f:
            binary = f.read()

        ctrl = Controller('template', ['asdf'])

        ctrl.file_input.file_selectors[0].value = binary
        ctrl.file_input._action_load()
        assert ctrl.file_input.norm_state == 'PpiA-FD'
        assert ctrl.file_input.norm_exposure == 0.0

        ctrl.file_input.norm_state = self.control[0]
        ctrl.file_input.norm_exposure = self.control[1]

        ctrl.file_input.exp_state = self.state
        assert ctrl.file_input.exp_exposures == [0.0, 0.167, 0.5, 1.0, 5.0, 10.0, 30.000002]
        ctrl.file_input._action_parse()

        assert isinstance(ctrl.series, KineticsSeries)
        assert isinstance(ctrl.peptides, PeptideMasterTable)

    def test_coverage(self):
        ctrl = Controller('template', ['asdf'])
        ctrl.series = self.series

    def test_initial_guesses(self):
        ctrl = Controller('template', ['asdf'])
        ctrl.series = self.series

        ctrl.fit_control._action_fit()
        assert 'half-life' in ctrl.sources.keys()