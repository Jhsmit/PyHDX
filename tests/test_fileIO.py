from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries
from pyhdx.fileIO import txt_to_protein, txt_to_np, read_dynamx, fmt_export, _get_f_width
from pyhdx.panel.apps import _main_app, _diff_app
from pyhdx.panel.data_sources import DataSource
from pathlib import Path
from io import StringIO
import pytest

directory = Path(__file__).parent

class TestFileIO(object):

    @classmethod
    def setup_class(cls):
        cls.fpath = directory / 'test_data' / 'ecSecB_apo.csv'


    def test_read_dynamx(self):
        data = read_dynamx(self.fpath)

        assert data.size == 567
        assert data['start'][0] == 9
        assert data['end'][0] == 18
        data = read_dynamx(self.fpath, intervals=('exclusive', 'inclusive'))
        assert data['start'][0] == 10

        data = read_dynamx(self.fpath, intervals=('inclusive', 'exclusive'))
        assert data['end'][0] == 17

        with pytest.raises(ValueError):
            data = read_dynamx(self.fpath, intervals=('foo', 'inclusive'))
        with pytest.raises(ValueError):
            data = read_dynamx(self.fpath, intervals=('inclusive', 'bar'))

        with open(self.fpath, mode='r') as f:
            read_dynamx(StringIO(f.read()))


    def test_methods(self):
        data = read_dynamx(self.fpath)
        path = directory / 'test_data' / 'ecSecB_guess.txt';

        #testing fmt_export
        with pytest.raises(ValueError):
            fmt_export(data, delimiter=',', width="foo")

        #testing txt_to_protein
        txt_to_protein(path)

        #testing txt_to_np
        with open(path, mode='r') as f:
            txt_to_np(StringIO('f'))


