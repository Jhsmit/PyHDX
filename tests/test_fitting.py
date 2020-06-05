import pytest
import os
from pyhdx import PeptideMeasurements, PeptideMasterTable, KineticsSeries
from pyhdx.fileIO import read_dynamx
from pyhdx.fitting import KineticsFitting, KineticsFitResult
from pyhdx.support import fmt_export, np_from_txt
import numpy as np

directory = os.path.dirname(__file__)
np.random.seed(43)


# class TestFitting(object):
#
#     def test_fit_section(self):
#         fpath = os.path.join(directory, 'test_data', 'ds1.csv')
#         control_100 = ('PpiA-FD', 0.167)
#         control_100 = ('PpiANative', 30.000002)
#
#         state = 'PpiANative'
#
#         data = read_dynamx(fpath)
#         print(np.unique(data['exposure']))
#         pcf = PeptideMasterTable(data)
#         #b = pcf.data['start'] < 50
#         #pcf.data = pcf.data[b]
#
#         pcf.set_control(control_100)
#         states = pcf.groupby_state()
#         series = states[state]
#
#         kf = KineticsFitting(series)
#         fr1 = kf.weighted_avg_fit()
#
#         arr1 = fr1.get_output(['rate', 'tau', 'tau1', 'tau2', 'r'])
#
#         gt = np_from_txt(os.path.join(directory, 'test_data', 'fit1_result.txt'))
#         for name in ['rate', 'tau', 'tau1', 'tau2', 'r']:
#             np.testing.assert_allclose(arr1[name], gt[name], rtol=1e-2)
#
#         fr2 = kf.global_fit(arr1)
#         arr2 = fr2.get_output(['rate', 'tau', 'tau1', 'tau2', 'r'])
#
#         gt = np_from_txt(os.path.join(directory, 'test_data', 'fit2_result.txt'))
#         for name in ['rate', 'tau', 'tau1', 'tau2', 'r']:
#             np.testing.assert_allclose(arr2[name], gt[name], rtol=1e-2)


class TestSimulatedDataFit(object):
    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'simulated_data.csv')
        cls.data = np_from_txt(fpath, delimiter=',')
        cls.sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

        cls.timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
        cls.start, cls.end = 5, 45 # total span of protein (inc, inc)
        cls.nc_start, cls.nc_end = 31, 34 # span of no coverage area (inc, inc)

    def test_fitting(self):
        pmt = PeptideMasterTable(self.data, drop_first=0, ignore_prolines=False, remove_nan=False)
        states = pmt.groupby_state()
        series = states['state1']

        kf = KineticsFitting(series, bounds=(1e-2, 800))

        fr1 = kf.weighted_avg_fit()
        out1 = fr1.get_output(['rate', 'k1', 'k2', 'r'])

        fr2 = kf.global_fit(out1)
        out2 = fr2.get_output(['rate', 'k1', 'k2', 'r'])

        check1 = np_from_txt(os.path.join(directory, 'test_data', 'Fit_simulated_wt_avg.txt'))
        check2 = np_from_txt(os.path.join(directory, 'test_data', 'Fit_simulated_global.txt'))

        for name in ['rate', 'k1', 'k2', 'r']:
            np.testing.assert_array_almost_equal(out1[name], check1[name], decimal=4)
            np.testing.assert_array_almost_equal(out2[name], check2[name], decimal=4)

