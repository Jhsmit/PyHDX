import pytest
import os
from pyhdx import PeptideMeasurements, PeptideCSVFile, KineticsSeries
import numpy as np

directory = os.path.dirname(__file__)


class TestUptakeFileModels(object):

    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        cls.pf = PeptideCSVFile(fpath)

    def test_peptidemeasurement(self):
        assert isinstance(self.pf, PeptideCSVFile)

        p_dict = self.pf.return_by_name('PpiA-FD', 0.167)

        name = 'PpiANative_0.167'
        pm = p_dict[name]
        assert pm.name == name
        assert isinstance(pm, PeptideMeasurements)

        res_scores = np.arange(pm.prot_len)
        scores = pm.calc_scores(res_scores)
        assert len(scores) == len(pm)

    def test_apply_controls(self):
        states_dict = self.pf.groupby_state()

        series1 = states_dict['PpiANative']
        assert len(series1) == 7
        assert len(series1.times) == 7
        assert isinstance(series1, KineticsSeries)
        control_100 = self.pf.get_data('PpiA-FD', 0.167)

        for pm in series1:
            assert len(pm) == 78
        series1.set_control(control_100)
        for pm in series1:
            assert len(pm) == 78

    def test_ds2(self):
        fpath = os.path.join(directory, 'test_data', 'ds2.csv')
        self.pf2 = PeptideCSVFile(fpath)

        states = self.pf2.groupby_state_control(('FD', 0.001), ('Native folded', 60.000004), remove_nan=False)
        series = states['folding_4C_10secLabelling']
        assert len(series[0]) == 80
        assert np.all(np.isnan(series[0].scores_average))

        states = self.pf2.groupby_state_control(('FD', 0.001), ('Native folded', 60.000004), remove_nan=True)
        series = states['folding_4C_10secLabelling']
        assert ~np.all(np.isnan(series[0].scores_average))
        assert len(series[0]) == 79
        
