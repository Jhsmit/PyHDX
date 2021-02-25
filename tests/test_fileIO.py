from pyhdx import PeptideMasterTable, read_dynamx, KineticsSeries
from pyhdx.fileIO import txt_to_protein, txt_to_np, read_dynamx
from pyhdx.panel.apps import _main_app, _diff_app
from pyhdx.panel.data_sources import DataSource
from pathlib import Path
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
            data = read_dynamx(self.fpath, intervals=('foo', 'bar'))
