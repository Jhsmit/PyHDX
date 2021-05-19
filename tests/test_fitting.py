import os
from pyhdx import PeptideMeasurements, PeptideMasterTable, KineticsSeries
from pyhdx.fileIO import read_dynamx, fmt_export, txt_to_protein, txt_to_np, csv_to_protein
from pyhdx.fitting import KineticsFitting, KineticsFitResult, BatchFitting, fit_rates_weighted_average, fit_gibbs_global
import numpy as np
import torch
import pandas as pd
import time
from dask.distributed import LocalCluster
import asyncio

directory = os.path.dirname(__file__)
np.random.seed(43)
torch.manual_seed(43)


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
        cls.series_apo = KineticsSeries(pf.get_state('SecB WT apo'), temperature=cls.temperature, pH=cls.pH)
        cls.series_dimer = KineticsSeries(pf.get_state('SecB his dimer apo'), temperature=cls.temperature, pH=cls.pH)

        data = pf.get_state('SecB WT apo')
        reduced_data = data[data['end'] < 40]
        cls.reduced_series = KineticsSeries(reduced_data)

        cluster = LocalCluster()
        cls.address = cluster.scheduler_address

    def test_initial_guess(self):
        result = fit_rates_weighted_average(self.reduced_series).compute()
        output = result.output

        assert output.size == 100

        # todo additional assertl, compare to stored values

    def test_global_fit(self):
        #kf = KineticsFitting(self.series_apo, bounds=(1e-2, 800), temperature=self.temperature, pH=self.pH)
        initial_rates = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))

        t0 = time.time()  # Very crude benchmarks
        gibbs_guess = self.series_apo.guess_deltaG(initial_rates['rate']).to_numpy()
        fr_global = fit_gibbs_global(self.series_apo, gibbs_guess, epochs=1000)
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

        output.to_file('compare.txt', fmt='pprint')

        check_protein = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_batch.csv'), column_depth=2)
        states = ['SecB WT apo', 'SecB his dimer apo']

        for state in states:
            from pandas.testing import assert_series_equal
            #
            result = output[state]['deltaG']
            test = check_protein[state]['deltaG']

            # # todo figure out the differences
            # idx = np.argmax(np.abs(result-test))
            # merged = pd.concat([result, test], axis=1)
            # merged['equal'] = (result - test).abs() > 0.01
            # merged['diff'] = (result - test).abs()
            #
            # from pyhdx.support import pprint_df_to_file
            # pprint_df_to_file(merged, 'merged.txt')

            assert_series_equal(result, test, rtol=0.1)


        result = asyncio.get_event_loop().run_until_complete(bf.global_fit_async(epochs=1000))
        output = result.output
        for state in states:
            result = output[state]['deltaG']
            test = check_protein[state]['deltaG']
            assert_series_equal(result, test, rtol=0.1)

            #assert np.allclose(output[state]['deltaG'], check_protein[state]['deltaG'], equal_nan=True, rtol=0.01)

    # def test_aligned_fit(self):

