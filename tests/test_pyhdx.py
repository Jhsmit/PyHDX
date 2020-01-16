import pytest
import os
from pyhdx import PeptideMeasurements, PeptideCSVFile

directory = os.path.dirname(__file__)


class TestUptakeFileModels(object):

    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        cls.pf = PeptideCSVFile(fpath)

    def test_peptidecsvfile(self):
        assert isinstance(self.pf, PeptideCSVFile)

        p_dict = self.pf.return_by_name('PpiA-FD', 0.167)

        pm = p_dict['PpiANative_0.167']
        assert isinstance(pm, PeptideMeasurements)
