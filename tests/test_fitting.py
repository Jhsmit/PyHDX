import pytest
import os
from pyhdx import PeptideMeasurements, PeptideCSVFile, KineticsSeries
from pyhdx.fileIO import read_dynamx
from pyhdx.fitting import KineticsFitting, KineticsFitResult
from pyhdx.support import fmt_export, np_from_txt
import numpy as np

directory = os.path.dirname(__file__)
np.random.seed(43)


class TestFitting(object):

    def test_fit_section(self):
        fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        control_100 = ('PpiA-FD', 0.167)
        control_100 = ('PpiANative', 30.000002)

        state = 'PpiANative'

        data = read_dynamx(fpath)
        print(np.unique(data['exposure']))
        pcf = PeptideCSVFile(data)
        #b = pcf.data['start'] < 50
        #pcf.data = pcf.data[b]

        pcf.set_control(control_100)
        states = pcf.groupby_state()
        series = states[state]

        kf = KineticsFitting(series)
        fr1 = kf.weighted_avg_fit()

        arr1 = fr1.get_output(['rate', 'tau', 'tau1', 'tau2', 'r'])

        gt = np_from_txt(os.path.join(directory, 'test_data', 'fit1_result.txt'))
        for name in ['rate', 'tau', 'tau1', 'tau2', 'r']:
            np.testing.assert_allclose(arr1[name], gt[name], rtol=1e-2)

        fr2 = kf.lsq_fit_blocks(arr1)
        arr2 = fr2.get_output(['rate', 'tau', 'tau1', 'tau2', 'r'])

        gt = np_from_txt(os.path.join(directory, 'test_data', 'fit2_result.txt'))
        for name in ['rate', 'tau', 'tau1', 'tau2', 'r']:
            np.testing.assert_allclose(arr2[name], gt[name], rtol=1e-2)


