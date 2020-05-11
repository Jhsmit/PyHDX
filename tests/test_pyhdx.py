import pytest
import os
from pyhdx import PeptideMeasurements, PeptideCSVFile, KineticsSeries
from pyhdx.fileIO import read_dynamx
import numpy as np
from functools import reduce
from operator import add

directory = os.path.dirname(__file__)


class TestUptakeFileModels(object):

    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        cls.pf1 = PeptideCSVFile(read_dynamx(fpath))

        fpath = os.path.join(directory, 'test_data', 'ds2.csv')
        cls.pf2 = PeptideCSVFile(read_dynamx(fpath))

        fpath = os.path.join(directory, 'test_data', 'ds3.csv')
        cls.pf3 = PeptideCSVFile(read_dynamx(fpath))

    def test_peptidemeasurement(self):
        assert isinstance(self.pf1, PeptideCSVFile)

        p_dict = self.pf1.return_by_name('PpiA-FD', 0.167)

        name = 'PpiANative_0.167'
        pm = p_dict[name]
        assert pm.name == name
        assert isinstance(pm, PeptideMeasurements)

        res_scores = np.arange(pm.prot_len)
        scores = pm.calc_scores(res_scores)
        assert len(scores) == len(pm)

    # def test_apply_controls(self):
    #     states_dict = self.pf1.groupby_state()
    #
    #     series1 = states_dict['PpiANative']
    #     assert len(series1) == 7
    #     assert len(series1.times) == 7
    #     assert isinstance(series1, KineticsSeries)
    #     control_100 = self.pf1.get_data('PpiA-FD', 0.167)
    #
    #     for pm in series1:
    #         assert len(pm) == 78
    #     series1.set_control(control_100)
    #     for pm in series1:
    #         assert len(pm) == 78

    def test_split(self):
        #control_100 = ("Full Deuteration control", 0.167)
        series_name = 'SecA-monomer'

        states = self.pf3.groupby_state()

        series = states[series_name]
        series.make_uniform()

        split_series = series.split()
        new_len = reduce(add, [reduce(add, [len(pm.data) for pm in ks]) for ks in split_series.values()])

        assert len(series.full_data) == new_len

    def test_uniform(self):
        control_100 = ("Full Deuteration control", 0.167)
        series_name = 'SecA-monomer'

        states = self.pf3.groupby_state()
        series = states[series_name]

        assert ~series.uniform
        # new_series = series.make_uniform(in_place=False)
        # assert new_series.uniform
        series.make_uniform()
        assert series.uniform

        l1 = len(series[0])
        series.make_uniform()
        l2 = len(series[0])
        assert l1 == l2

    def test_ds2(self):
        self.pf2.set_control(('FD', 0.001), ('Native folded', 60.000004), remove_nan=False)
        states = self.pf2.groupby_state()
        series = states['folding_4C_10secLabelling']
        assert len(series[0]) == 80
        assert np.all(np.isnan(series[0].scores_average))

        self.pf2.set_control(('FD', 0.001), ('Native folded', 60.000004), remove_nan=True)
        states = self.pf2.groupby_state()
        series = states['folding_4C_10secLabelling']
        assert ~np.all(np.isnan(series[0].scores_average))
        assert len(series[0]) == 79

