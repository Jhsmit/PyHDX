import os
from pyhdx import PeptideMeasurements, PeptideMasterTable, KineticsSeries
from pyhdx.fileIO import read_dynamx, fmt_export, txt_to_protein, txt_to_np, csv_to_protein
from pyhdx.fitting import KineticsFitting, KineticsFitResult, BatchFitting
import numpy as np
import torch
import pandas as pd
import time
from dask.distributed import LocalCluster
import asyncio

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
            assert np.allclose(out1[name], check1[name], rtol=0.01, equal_nan=True)

    def test_torch_fitting(self):
        pmt = PeptideMasterTable(self.data, drop_first=1, ignore_prolines=True, remove_nan=False)
        pmt.set_backexchange(0.)
        states = pmt.groupby_state()
        series = states['state1']

        kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=300, pH=8)
        initial_rates = txt_to_protein(os.path.join(directory, 'test_data', 'fit_simulated_wt_avg.txt'))

        fr_pfact = kf.global_fit(initial_rates, epochs=1000)
        out_deltaG = fr_pfact.output
        check_deltaG = txt_to_protein(os.path.join(directory, 'test_data', 'fit_simulated_torch.txt'))

        assert np.allclose(check_deltaG['deltaG'], out_deltaG['deltaG'], equal_nan=True, rtol=0.01)


class TestSecBDataFit(object):
    @classmethod
    def setup_class(cls):
        fpath_apo = os.path.join(directory, 'test_data', 'ecSecB_apo.csv')
        fpath_dimer = os.path.join(directory, 'test_data', 'ecSecB_dimer.csv')
        data = read_dynamx(fpath_apo, fpath_dimer)
        control = ('Full deuteration control', 0.167)

        cls.temperature, cls.pH = 273.15 + 30, 8.

        pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
        pf.set_control(control)
        states = pf.groupby_state()
        cls.series_apo = states['SecB WT apo']
        cls.series_dimer = states['SecB his dimer apo']

        cluster = LocalCluster()
        cls.address = cluster.scheduler_address

    def test_global_fit(self):
        kf = KineticsFitting(self.series_apo, bounds=(1e-2, 800), temperature=self.temperature, pH=self.pH)
        initial_rates = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))

        t0 = time.time()  # Very crude benchmarks
        fr_global = kf.global_fit(initial_rates, epochs=1000)
        t1 = time.time()

        assert t1 - t0 < 5
        out_deltaG = fr_global.output
        check_deltaG = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_torch_fit.txt'))

        assert np.allclose(check_deltaG['deltaG'], out_deltaG['deltaG'], equal_nan=True, rtol=0.01)
        assert np.allclose(check_deltaG['covariance'], out_deltaG['covariance'], equal_nan=True, rtol=0.01)

    def test_global_fit_async(self):
        kf = KineticsFitting(self.series_apo, bounds=(1e-2, 800), temperature=self.temperature, pH=self.pH, cluster=self.address)
        initial_rates = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))

        t0 = time.time()  # Very crude benchmarks
        fr_global = asyncio.get_event_loop().run_until_complete(kf.global_fit_async(initial_rates, epochs=1000))  # py37: run
        t1 = time.time()

        out_deltaG = fr_global.output
        check_deltaG = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_torch_fit.txt'))

        assert np.allclose(check_deltaG['deltaG'], out_deltaG['deltaG'], equal_nan=True, rtol=0.01)
        assert np.allclose(check_deltaG['covariance'], out_deltaG['covariance'], equal_nan=True, rtol=0.01)

    def test_batch_fit(self):
        kfs = [KineticsFitting(series, temperature=self.temperature, pH=self.pH) for series in [self.series_apo, self.series_dimer]]
        guess = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))

        bf = BatchFitting(kfs, [guess, guess])
        result = bf.global_fit(epochs=1000)
        output = result.output

        check_protein = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_batch.csv'))
        states = ['SecB WT apo', 'SecB his dimer apo']

        for state in states:
            assert np.allclose(output[state]['deltaG'], check_protein[state]['deltaG'], equal_nan=True, rtol=0.01)

        result = asyncio.get_event_loop().run_until_complete(bf.global_fit_async(epochs=1000))
        output = result.output
        for state in states:
            assert np.allclose(output[state]['deltaG'], check_protein[state]['deltaG'], equal_nan=True, rtol=0.01)

