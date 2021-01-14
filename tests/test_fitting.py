import os
from pyhdx import PeptideMeasurements, PeptideMasterTable, KineticsSeries
from pyhdx.fileIO import read_dynamx, fmt_export, txt_to_protein, txt_to_np
from pyhdx.fitting import KineticsFitting, KineticsFitResult
import numpy as np
import torch

directory = os.path.dirname(__file__)
np.random.seed(43)
torch.manual_seed(43)

class TestSimulatedDataFit(object):
    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'simulated_data_uptake.csv')
        cls.data = txt_to_np(fpath, delimiter=',')
        cls.data['end'] += 1  # because this simulated data is in old format of inclusive, inclusive
        cls.sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

        cls.timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
        cls.start, cls.end = 5, 46 # total span of protein (inc, ex)
        cls.nc_start, cls.nc_end = 31, 35 # span of no coverage area (inc, ex)

    def test_fitting(self):
        pmt = PeptideMasterTable(self.data, drop_first=1, ignore_prolines=True, remove_nan=False)
        pmt.set_backexchange(0.)
        states = pmt.groupby_state()
        series = states['state1']

        kf = KineticsFitting(series, bounds=(1e-2, 800))
        fr1 = kf.weighted_avg_fit()

        out1 = fr1.output
        check1 = txt_to_protein(os.path.join(directory, 'test_data', 'fit_simulated_wt_avg.txt'))
        for name in ['rate', 'k1', 'k2', 'r']:
            np.allclose(out1[name], check1[name])

    def test_torch_fitting(self):
        pmt = PeptideMasterTable(self.data, drop_first=1, ignore_prolines=True, remove_nan=False)
        pmt.set_backexchange(0.)
        states = pmt.groupby_state()
        series = states['state1']

        kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=300, pH=8)
        initial_rates = txt_to_protein(os.path.join(directory, 'test_data', 'fit_simulated_wt_avg.txt'))

        fr_pfact = kf.global_fit_torch(initial_rates, epochs=1000)
        out_deltaG = fr_pfact.output
        check_deltaG = txt_to_protein(os.path.join(directory, 'test_data', 'fit_simulated_torch.txt'))

        np.allclose(check_deltaG['deltaG'], out_deltaG['deltaG'], equal_nan=True)

