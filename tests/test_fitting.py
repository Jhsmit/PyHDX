import os
from pyhdx import PeptideMasterTable, HDXMeasurement
from pyhdx.fileIO import read_dynamx, csv_to_protein
from pyhdx.fitting import fit_rates_weighted_average, fit_gibbs_global, fit_gibbs_global_batch, fit_gibbs_global_batch_aligned
from pyhdx.models import HDXMeasurementSet
import numpy as np
import torch
import time
from dask.distributed import LocalCluster

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
        cls.series_apo = HDXMeasurement(pf.get_state('SecB WT apo'), temperature=cls.temperature, pH=cls.pH)
        cls.series_dimer = HDXMeasurement(pf.get_state('SecB his dimer apo'), temperature=cls.temperature, pH=cls.pH)

        data = pf.get_state('SecB WT apo')
        reduced_data = data[data['end'] < 40]
        cls.reduced_series = HDXMeasurement(reduced_data)

        cluster = LocalCluster()
        cls.address = cluster.scheduler_address

    def test_initial_guess(self):
        result = fit_rates_weighted_average(self.reduced_series)
        output = result.output

        assert output.size == 100

        # todo additional tests:
        #  result = fit_rates_half_time_interpolate()

        # todo additional assert, compare to stored values

    def test_global_fit(self):
        #kf = KineticsFitting(self.series_apo, bounds=(1e-2, 800), temperature=self.temperature, pH=self.pH)
        initial_rates = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))

        t0 = time.time()  # Very crude benchmarks
        gibbs_guess = self.series_apo.guess_deltaG(initial_rates['rate']).to_numpy()
        fr_global = fit_gibbs_global(self.series_apo, gibbs_guess, epochs=1000, r1=2)
        t1 = time.time()

        assert t1 - t0 < 5
        out_deltaG = fr_global.output
        check_deltaG = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_torch_fit.txt'))

        assert np.allclose(check_deltaG['deltaG'], out_deltaG['deltaG'], equal_nan=True, rtol=0.01)
        assert np.allclose(check_deltaG['covariance'], out_deltaG['covariance'], equal_nan=True, rtol=0.01)

    def test_batch_fit(self):
        hdx_set = HDXMeasurementSet([self.series_apo, self.series_dimer])
        guess = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_guess.txt'))

        gibbs_guess = hdx_set.guess_deltaG([guess['rate'], guess['rate']])
        result = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=1000)

        output = result.output

        check_protein = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_batch.csv'), column_depth=2)
        states = ['SecB WT apo', 'SecB his dimer apo']

        for state in states:
            from pandas.testing import assert_series_equal

            result = output[state]['deltaG']
            test = check_protein[state]['deltaG']

            assert_series_equal(result, test, rtol=0.1)

        mock_alignment = {
            'apo': 'MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA',
            'dimer': 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------',
        }

        hdx_set.add_alignment(list(mock_alignment.values()))

        gibbs_guess = hdx_set.guess_deltaG([guess['rate'], guess['rate']])
        aligned_result = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000)
        output = aligned_result.output
        check_protein = csv_to_protein(os.path.join(directory, 'test_data', 'ecSecB_batch_aligned.csv'), column_depth=2)
        states = ['SecB WT apo', 'SecB his dimer apo']

        for state in states:
            from pandas.testing import assert_series_equal
            result = output[state]['deltaG']
            test = check_protein[state]['deltaG']

            assert_series_equal(result, test, rtol=0.1)


