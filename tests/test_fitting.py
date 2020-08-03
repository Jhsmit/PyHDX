import pytest
import os
from pyhdx import PeptideMeasurements, PeptideMasterTable, KineticsSeries
from pyhdx.fileIO import read_dynamx
from pyhdx.fitting import KineticsFitting, KineticsFitResult
from pyhdx.support import fmt_export, np_from_txt
import numpy as np

import matplotlib.pyplot as plt

directory = os.path.dirname(__file__)
np.random.seed(43)


class TestDissociationFitting(object):

    def test_fit_section(self):
        pass
        # fpath = os.path.join(directory, 'test_data', 'ds2.csv')
        # control_100 = ('FD', 0.001)
        # control_0 = ('Native folded', 60.000004)
        #
        # state = 'folding_4C_10secLabelling'
        # data = read_dynamx(fpath)
        # pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True)
        # pf.set_control(control_100, control_0)
        # states = pf.groupby_state()
        # series = states[state]
        # split = list(series.split().items())[-1]
        #
        # kf = KineticsFitting(split)
        # fr1 = kf.weighted_avg_fit(model_type='dissociation')
        # arr1 = fr1.get_output(['rate', 'k1', 'k2', 'r'])
        #
        # fr2 = kf.blocks_fit(arr1, model_type='dissociation')
        # arr2 = fr2.get_output(['rate', 'k1', 'k2', 'r'])


class TestSimulatedDataFit(object):
    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'simulated_data.csv')
        cls.data = np_from_txt(fpath, delimiter=',')
        cls.data['end'] += 1  # because this simulated data is in old format of inclusive, inclusive
        cls.sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

        cls.timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
        cls.start, cls.end = 5, 46 # total span of protein (inc, ex)
        cls.nc_start, cls.nc_end = 31, 35 # span of no coverage area (inc, ex)

    def test_fitting(self):
        np.random.seed(43)
        pmt = PeptideMasterTable(self.data, drop_first=1, ignore_prolines=True, remove_nan=False)
        states = pmt.groupby_state()
        series = states['state1']

        kf = KineticsFitting(series, bounds=(1e-2, 800))

        fr1 = kf.weighted_avg_fit()
        out1 = fr1.get_output(['rate', 'k1', 'k2', 'r'])

        fr2 = kf.blocks_fit(out1)
        out2 = fr2.get_output(['rate', 'k1', 'k2', 'r'])

        check1 = np_from_txt(os.path.join(directory, 'test_data', 'Fit_simulated_wt_avg.txt'))
        check2 = np_from_txt(os.path.join(directory, 'test_data', 'Fit_simulated_blocks.txt'))

        for name in ['rate', 'k1', 'k2', 'r']:
            #indices = np.searchsorted(out1['r_number'], check1['r_number'])

            np.testing.assert_array_almost_equal(out1[name], check1[name], decimal=4)
            np.testing.assert_array_almost_equal(out2[name], check2[name], decimal=4)

    def test_tf_pfact_fitting(self):
        pmt = PeptideMasterTable(self.data, drop_first=1, ignore_prolines=True, remove_nan=False)
        states = pmt.groupby_state()
        series = states['state1']

        kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=300, pH=8)
        initial_rates = np_from_txt(os.path.join(directory, 'test_data', 'Fit_simulated_wt_avg.txt'))

        fr_pfact = kf.global_fit_new(initial_rates, use_kint=True)
        out_pfact = fr_pfact.output

        check_pfact = np_from_txt(os.path.join(directory, 'test_data', 'Fit_simulated_pfact.txt'))
        indices = np.searchsorted(out_pfact['r_number'], check_pfact['r_number'])

        np.testing.assert_array_almost_equal(out_pfact['log_P'][indices], check_pfact['log_P'], decimal=1)
