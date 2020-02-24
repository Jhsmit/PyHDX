import pytest
import os
from pyhdx import PeptideMeasurements, PeptideCSVFile, KineticsSeries
from pyhdx.fitting import KineticsFitting, KineticsFitResult
import numpy as np

directory = os.path.dirname(__file__)


class TestFitting(object):

    # @classmethod
    # def setup_class(cls):


    def test_fit_section(self):
        fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        control_100 = ('PpiA-FD', 0.167)
        state = 'PpiANative'

        pcf = PeptideCSVFile(fpath)
        b = pcf.data['start'] < 50
        pcf.data = pcf.data[b]

        states = pcf.groupby_state_control(control_100)
        series = states[state]

        kf = KineticsFitting(series)
        r, m, b = kf.do_fitting()
        fr = KineticsFitResult(r, m, b)

        rate = fr.rate
        expected = np.array([1.00335533, 1.00335533, 1.00335533, 1.00335533, 1.00335533, 1.00335533,
                             1.00335533, 1.00335533, 0.65508649, 0.63973347, 0.70839084, 0.70839084,
                             0.70839084, 0.70839084, 0.70839084, 0.66683659, 0.60067197, 0.43035482,
                             0.33043168, 0.30382436, 0.29411784, 0.27600381, 0.27600381, 0.27600381,
                             0.27600381, 0.27600381, 0.27600381, 0.27600381, 0.27600381, 0.27600381,
                             0.27600381, 0.27418457, 0.26622289, 0.26622289, 0.26759417, 0.27282912,
                             0.27282912, 0.27282912, 0.27282912, 0.27282912, 0.23271286, 0.23271286,
                             0.23271286, 0.23271286, 0.23271286, 0.23271286, 0.21722117, 0.21722117,
                             0.21722117, 0.21722117, 0.21722117, 0.21722117, 0.21722117, 0.21722117])

        np.testing.assert_allclose(expected, rate, rtol=1e-7)
