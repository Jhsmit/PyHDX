from pyhdx.fileIO import csv_to_protein, txt_to_np, read_dynamx, fmt_export, csv_to_np
from pyhdx.models import Protein
from pathlib import Path
from io import StringIO
import numpy as np
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
            data = read_dynamx(StringIO(f.read()))
            assert data.size == 567

    def test_fmt_export(self):
        # testing fmt_export
        data = read_dynamx(self.fpath)
        with pytest.raises(ValueError):
            fmt_export(data, delimiter=',', width="foo")

        fmt, hdr = fmt_export(data, delimiter=',', width=0)
        assert 'protein' in hdr

        fmt, hdr = fmt_export(data, delimiter=',', width='auto')
        assert 'protein' in hdr

        fmt, hdr = fmt_export(data, header=False)
        assert len(hdr) == 0

        data = np.genfromtxt(self.fpath, skip_header=1, delimiter=',', dtype=[('foo', 'V'), ('bar', 'V')])
        with pytest.raises(TypeError):
            fmt_export(data, delimiter=',', width=0)

    # testing other functions in fileIO
    def test_methods(self):
        path = directory / 'test_data' / 'ecSecB_guess.txt'

        # testing csv_to_protein
        ret = csv_to_protein(path)
        assert type(ret) == Protein
        assert ret.index.name == 'r_number'

        # testing txt_to_np
        with open(path, mode='r') as f:
            ret = csv_to_np(StringIO(f.read()))
            assert 'r_number' in ret.dtype.names
