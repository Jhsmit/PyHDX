import os
import tempfile

from pyhdx import PeptideMasterTable, HDXMeasurement
from pyhdx.fileIO import read_dynamx, csv_to_protein, csv_to_dataframe, save_fitresult, load_fitresult
from pyhdx.fitting import fit_rates_weighted_average, fit_gibbs_global, fit_gibbs_global_batch, fit_gibbs_global_batch_aligned
from pyhdx.models import HDXMeasurementSet
import numpy as np
import torch
import time
from dask.distributed import LocalCluster
from pathlib import Path

import pandas as pd

directory = Path(__file__).parent
np.random.seed(43)
torch.manual_seed(43)


class TestSecBDataFit(object):
    @classmethod
    def setup_class(cls):
        fpath_apo = directory / 'test_data' / 'ecSecB_apo.csv'
        fpath_dimer = directory / 'test_data' / 'ecSecB_dimer.csv'
        data = read_dynamx(fpath_apo, fpath_dimer)
        control = ('Full deuteration control', 0.167)

        cls.temperature, cls.pH = 273.15 + 30, 8.

        pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
        pf.set_control(control)
        cls.hdxm_apo = HDXMeasurement(pf.get_state('SecB WT apo'), temperature=cls.temperature, pH=cls.pH)
        cls.hdxm_dimer = HDXMeasurement(pf.get_state('SecB his dimer apo'), temperature=cls.temperature, pH=cls.pH)

        data = pf.get_state('SecB WT apo')
        reduced_data = data[data['end'] < 40]
        cls.reduced_hdxm = HDXMeasurement(reduced_data)

        cluster = LocalCluster()
        cls.address = cluster.scheduler_address

    def test_initial_guess(self):
        result = fit_rates_weighted_average(self.reduced_hdxm)
        output = result.output

        assert output.size == 100
        check_rates = csv_to_protein(directory / 'test_data' / 'ecSecB_reduced_guess.csv')
        pd.testing.assert_series_equal(check_rates['rate'], output['rate'])

        # todo additional tests:
        #  result = fit_rates_half_time_interpolate()

    def test_global_fit(self):
        # initial_rates = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))
        initial_rates = pd.read_csv(directory / 'test_data' / 'ecSecB_guess.txt', header=[0], comment='#', index_col=0)

        t0 = time.time()  # Very crude benchmarks
        gibbs_guess = self.hdxm_apo.guess_deltaG(initial_rates['rate']).to_numpy()
        fr_global = fit_gibbs_global(self.hdxm_apo, gibbs_guess, epochs=1000, r1=2)
        t1 = time.time()

        assert t1 - t0 < 5
        out_deltaG = fr_global.output
        check_deltaG = csv_to_protein(directory / 'test_data' / 'ecSecB_torch_fit.csv')

        assert np.allclose(check_deltaG['deltaG'], out_deltaG['deltaG'], equal_nan=True, rtol=0.01)
        assert np.allclose(check_deltaG['covariance'], out_deltaG['covariance'], equal_nan=True, rtol=0.01)
        assert np.allclose(check_deltaG['k_obs'], out_deltaG['k_obs'], equal_nan=True, rtol=0.01)

        #These should perhaps be moved to fileIO tests
        with tempfile.TemporaryDirectory() as tempdir:
            fpath = Path(tempdir) / 'fit_result_single.csv'
            fr_global.to_file(fpath)
            df = csv_to_dataframe(fpath)
            assert df.attrs['metadata'] == fr_global.metadata

            fit_result_dir = Path(tempdir) / 'fit_result'
            save_fitresult(fit_result_dir, fr_global, log_lines=['test123'])

            log_lines = Path(fit_result_dir / 'log.txt').read_text().split('\n')
            assert log_lines[-1] == 'test123'

            fit_result_loaded = load_fitresult(fit_result_dir)
            assert isinstance(fit_result_loaded.losses, pd.DataFrame)
            assert isinstance(fit_result_loaded.data_obj, HDXMeasurement)

            timepoints = np.linspace(0, 30, num=100)
            d_calc = fit_result_loaded(timepoints)


    def test_batch_fit(self):
        hdx_set = HDXMeasurementSet([self.hdxm_apo, self.hdxm_dimer])
        #guess = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))
        guess = pd.read_csv(directory / 'test_data' / 'ecSecB_guess.txt', header=[0], comment='#', index_col=0)

        gibbs_guess = hdx_set.guess_deltaG([guess['rate'], guess['rate']])
        result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=1000)

        with tempfile.TemporaryDirectory() as tempdir:
            fpath = Path(tempdir) / 'fit_result_batch.csv'
            result.to_file(fpath)
            df = csv_to_dataframe(fpath)
            assert df.attrs['metadata'] == result.metadata

        output = result.output

        check_protein = csv_to_protein(directory / 'test_data' / 'ecSecB_batch.csv')
        states = ['SecB WT apo', 'SecB his dimer apo']

        for state in states:
            from pandas.testing import assert_series_equal

            result = output[state]['deltaG']
            test = check_protein[state]['deltaG']

            assert_series_equal(result, test, rtol=0.1)

        mock_alignment = {
            'apo':   'MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA',
            'dimer': 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------',
        }

        hdx_set.add_alignment(list(mock_alignment.values()))

        gibbs_guess = hdx_set.guess_deltaG([guess['rate'], guess['rate']])
        aligned_result = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000)
        output = aligned_result.output
        check_protein = csv_to_protein(directory / 'test_data' / 'ecSecB_batch_aligned.csv')
        states = ['SecB WT apo', 'SecB his dimer apo']

        for state in states:
            from pandas.testing import assert_series_equal
            result = output[state]['deltaG']
            test = check_protein[state]['deltaG']

            assert_series_equal(result, test, rtol=0.1)
