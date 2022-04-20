import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch
import yaml
from pandas.testing import assert_series_equal, assert_frame_equal
from pyhdx import PeptideMasterTable, HDXMeasurement
from pyhdx.config import cfg
from pyhdx.fileIO import read_dynamx, csv_to_protein, csv_to_dataframe
from pyhdx.fitting import (
    fit_rates_weighted_average,
    fit_gibbs_global,
    fit_gibbs_global_batch,
    fit_gibbs_global_batch_aligned,
    fit_rates_half_time_interpolate,
    GenericFitResult,
)
from pyhdx.batch_processing import StateParser
from pyhdx.models import HDXMeasurementSet

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"

np.random.seed(43)
torch.manual_seed(43)

sequence = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"
sequence_dimer = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"


class TestSecBDataFit(object):
    @classmethod
    def setup_class(cls):
        fpath_apo = input_dir / "ecSecB_apo.csv"
        fpath_dimer = input_dir / "ecSecB_dimer.csv"
        data = read_dynamx(fpath_apo, fpath_dimer)
        control = ("Full deuteration control", 0.167 * 60)

        cls.temperature, cls.pH = 273.15 + 30, 8.0

        pf = PeptideMasterTable(
            data, drop_first=1, ignore_prolines=True, remove_nan=False
        )
        pf.set_control(control)
        cls.hdxm_apo = HDXMeasurement(
            pf.get_state("SecB WT apo"),
            temperature=cls.temperature,
            pH=cls.pH,
            sequence=sequence,
        )
        cls.hdxm_dimer = HDXMeasurement(
            pf.get_state("SecB his dimer apo"),
            temperature=cls.temperature,
            pH=cls.pH,
            sequence=sequence_dimer,
        )

        data = pf.get_state("SecB WT apo")
        reduced_data = data[data["end"] < 40]
        cls.reduced_hdxm = HDXMeasurement(reduced_data)

    def test_initial_guess_wt_average(self):
        result = fit_rates_weighted_average(self.reduced_hdxm)
        output = result.output

        assert output.size == 100
        check_rates = csv_to_protein(output_dir / "ecSecB_reduced_guess.csv")
        pd.testing.assert_series_equal(check_rates["rate"], output["rate"])

    def test_initial_guess_half_time_interpolate(self):
        result = fit_rates_half_time_interpolate(self.reduced_hdxm)
        assert isinstance(result, GenericFitResult)
        assert result.output.index.name == "r_number"
        assert result.output["rate"].mean() == pytest.approx(0.04343354509254464)

        # todo additional tests:
        #  result = fit_rates_half_time_interpolate()

    @pytest.mark.skip(reason="Hangs on GitHub Actions on Ubuntu")
    def test_dtype_cuda(self):
        check_deltaG = csv_to_protein(output_dir / "ecSecB_torch_fit.csv")
        initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

        cfg.set("fitting", "device", "cuda")
        gibbs_guess = self.hdxm_apo.guess_deltaG(initial_rates["rate"]).to_numpy()

        if torch.cuda.is_available():
            fr_global = fit_gibbs_global(self.hdxm_apo, gibbs_guess, epochs=1000, r1=2)
            out_deltaG = fr_global.output
            for field in ["dG", "k_obs", "covariance"]:
                assert_series_equal(
                    check_deltaG[field],
                    out_deltaG[self.hdxm_apo.name, field],
                    rtol=0.01,
                    check_dtype=False,
                    check_names=False,
                )
        else:
            with pytest.raises(AssertionError, match=r".* CUDA .*"):
                fr_global = fit_gibbs_global(
                    self.hdxm_apo, gibbs_guess, epochs=1000, r1=2
                )

        cfg.set("fitting", "device", "cpu")
        cfg.set("fitting", "dtype", "float32")

        fr_global = fit_gibbs_global(self.hdxm_apo, gibbs_guess, epochs=1000, r1=2)
        dg = fr_global.model.dG
        assert dg.dtype == torch.float32

        out_deltaG = fr_global.output
        for field in ["dG", "k_obs"]:
            assert_series_equal(
                check_deltaG[field],
                out_deltaG[self.hdxm_apo.name, field],
                rtol=0.01,
                check_dtype=False,
                check_names=False,
            )

        cfg.set("fitting", "dtype", "float64")

    def test_global_fit(self):
        initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

        t0 = time.time()  # Very crude benchmarks
        gibbs_guess = self.hdxm_apo.guess_deltaG(initial_rates["rate"])
        fr_global = fit_gibbs_global(self.hdxm_apo, gibbs_guess, epochs=1000, r1=2)
        t1 = time.time()

        # assert t1 - t0 < 5  # Fails sometimes
        out_deltaG = fr_global.output
        check_deltaG = csv_to_protein(output_dir / "ecSecB_torch_fit.csv")

        for field in ["dG", "covariance", "k_obs"]:
            assert_series_equal(
                check_deltaG[self.hdxm_apo.name, field],
                out_deltaG[self.hdxm_apo.name, field],
                rtol=0.01,
                check_names=False,
            )

        errors = fr_global.get_squared_errors()
        assert errors.shape == (1, self.hdxm_apo.Np, self.hdxm_apo.Nt)

    @pytest.mark.skip(
        reason="Longer fit is not checked by default due to long computation times"
    )
    def test_global_fit_extended(self):
        check_deltaG = csv_to_protein(output_dir / "ecSecB_torch_fit_epochs_20000.csv")
        initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")
        gibbs_guess = self.hdxm_apo.guess_deltaG(initial_rates["rate"])

        t0 = time.time()  # Very crude benchmarks
        fr_global = fit_gibbs_global(self.hdxm_apo, gibbs_guess, epochs=20000, r1=2)
        t1 = time.time()

        assert t1 - t0 < 50
        out_deltaG = fr_global.output
        for field in ["dG", "k_obs", "covariance"]:
            assert_series_equal(
                check_deltaG[self.hdxm_apo.name, field],
                out_deltaG[self.hdxm_apo.name, field],
                rtol=0.01,
                check_dtype=False,
            )

        errors = fr_global.get_squared_errors()
        assert errors.shape == (1, self.hdxm_apo.Np, self.hdxm_apo.Nt)

    @pytest.mark.skip(
        reason="Longer fit is not checked by default due to long computation times"
    )
    def test_global_fit_extended_cuda(self):
        check_deltaG = csv_to_protein(output_dir / "ecSecB_torch_fit_epochs_20000.csv")
        initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")
        gibbs_guess = self.hdxm_apo.guess_deltaG(initial_rates["rate"])

        # todo allow contextmanger?
        cfg.set("fitting", "device", "cuda")
        cfg.set("fitting", "dtype", "float32")

        fr_global = fit_gibbs_global(self.hdxm_apo, gibbs_guess, epochs=20000, r1=2)
        out_deltaG = fr_global.output

        for field in ["dG", "k_obs"]:
            assert_series_equal(
                check_deltaG[self.hdxm_apo.name, field],
                out_deltaG[self.hdxm_apo.name, field],
                rtol=0.01,
                check_dtype=False,
            )

        cfg.set("fitting", "device", "cpu")
        cfg.set("fitting", "dtype", "float64")

    def test_batch_fit(self, tmp_path):
        hdx_set = HDXMeasurementSet([self.hdxm_dimer, self.hdxm_apo])
        guess = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

        # Create rates dataframe
        rates_df = pd.DataFrame({name: guess["rate"] for name in hdx_set.names})
        gibbs_guess = hdx_set.guess_deltaG(rates_df)
        fr_global = fit_gibbs_global_batch(hdx_set, gibbs_guess, epochs=1000)

        fpath = Path(tmp_path) / "fit_result_batch.csv"
        fr_global.to_file(fpath)
        df = csv_to_dataframe(fpath)
        assert df.attrs["metadata"] == fr_global.metadata

        output = fr_global.output

        check_protein = csv_to_protein(output_dir / "ecSecB_batch.csv")
        states = ["SecB WT apo", "SecB his dimer apo"]

        for state in states:
            from pandas.testing import assert_series_equal

            result = output[state]["dG"]
            test = check_protein[state]["dG"]

            assert_series_equal(result, test, rtol=0.1)

        errors = fr_global.get_squared_errors()
        assert errors.shape == (hdx_set.Ns, hdx_set.Np, hdx_set.Nt)

        test = fr_global.get_peptide_mse().fillna(-1)
        ref = csv_to_dataframe(output_dir / "ecSecB_batch_peptide_mse.csv").fillna(-1)
        assert_frame_equal(test, ref, atol=1e-1, rtol=5e-1)

        test = fr_global.get_residue_mse().fillna(-1)
        ref = csv_to_dataframe(output_dir / "ecSecB_batch_residue_mse.csv").fillna(-1)
        assert_frame_equal(test, ref, atol=1e-1, rtol=5e-1)

        test = fr_global.losses.fillna(-1)
        ref = csv_to_dataframe(output_dir / "ecSecB_batch_loss.csv").fillna(-1)
        assert_frame_equal(test, ref, atol=1e-3, rtol=1e-2)

        # test alignment fit
        mock_alignment = {
            "apo": "MSEQNNTEMTFQIQRIYTKDI------------SFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG-------------------EETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRG----TFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA",
            "dimer": "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVY--------------EVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGA----YCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA-----------------",
        }

        hdx_set.add_alignment(list(mock_alignment.values()))

        gibbs_guess = hdx_set[0].guess_deltaG(
            guess["rate"]
        )  # Guesses from first measurement
        aligned_result = fit_gibbs_global_batch_aligned(
            hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000
        )
        output = aligned_result.output
        check_protein = csv_to_protein(output_dir / "ecSecB_batch_aligned.csv")
        states = ["SecB WT apo", "SecB his dimer apo"]

        for state in states:
            from pandas.testing import assert_series_equal

            result = output[state]["dG"]
            test = check_protein[state]["dG"]

            assert_series_equal(result, test, rtol=0.1)

    # batch fit on delta N/C tail dataset
    def test_batch_fit_delta(self, tmp_path):
        yaml_file = input_dir / "data_states_deltas.yaml"
        yaml_spec = yaml.safe_load(yaml_file.read_text())
        parser = StateParser(yaml_spec, data_src=input_dir)

        hdxm_set = parser.load_hdxmset()
        guess_output = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

        gibbs_guess = hdxm_set[0].guess_deltaG(guess_output["rate"])

        # broadcast single guess over samples
        fr_global = fit_gibbs_global_batch(hdxm_set, gibbs_guess, epochs=200)
        output = fr_global.output

        check = csv_to_dataframe(output_dir / "ecsecb_delta_batch" / "fit_result.csv")
        states = check.columns.unique(level=0)

        for state in states:
            from pandas.testing import assert_series_equal

            result = output[state]["dG"]
            test = check[state]["dG"]

            assert_series_equal(result, test, rtol=0.1)

        errors = fr_global.get_squared_errors()
        assert errors.shape == (hdxm_set.Ns, hdxm_set.Np, hdxm_set.Nt)
        assert not np.any(np.isnan(errors))
