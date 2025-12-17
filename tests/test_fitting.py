import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch
from hdxms_datasets import HDXDataSet, load_dataset
from pandas.testing import assert_frame_equal, assert_series_equal

from pyhdx import HDXMeasurement
from pyhdx.config import cfg
from pyhdx.fileIO import csv_to_dataframe, read_dynamx
from pyhdx.fitting import (
    GenericFitResult,
    fit_d_uptake,
    fit_gibbs_global,
    fit_gibbs_global_batch,
    fit_gibbs_global_batch_aligned,
    fit_rates_half_time_interpolate,
    fit_rates_weighted_average,
)
from pyhdx.models import HDXMeasurementSet
from pyhdx.process import apply_control, correct_d_uptake
from pyhdx.legacy import filter_peptides


cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"

np.random.seed(43)
torch.manual_seed(43)

sequence = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"
sequence_dimer = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPAARECIASMVARGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"

TEMPERATURE, PH = 273.15 + 30, 8.0


@pytest.fixture()
def dataframes() -> tuple[pd.DataFrame, pd.DataFrame]:
    df_apo = read_dynamx(input_dir / "ecSecB_apo.csv")
    df_dimer = read_dynamx(input_dir / "ecSecB_dimer.csv")

    fd = {"state": "Full deuteration control", "exposure": {"value": 0.167, "unit": "min"}}
    fd_df = filter_peptides(df_apo, **fd)

    apo_peptides = filter_peptides(df_apo, state="SecB WT apo")
    dimer_peptides = filter_peptides(df_dimer, state="SecB his dimer apo")

    apo_control = apply_control(apo_peptides, fd_df)
    dimer_control = apply_control(dimer_peptides, fd_df)

    apo_corrected = correct_d_uptake(apo_control)
    dimer_corrected = correct_d_uptake(dimer_control)

    return apo_corrected, dimer_corrected


@pytest.fixture()
def base_dataset() -> HDXDataSet:
    dataset_dir = input_dir / "HDX_D9096080"
    dataset = load_dataset(dataset_dir)

    return dataset


@pytest.fixture()
def hdxm_apo(dataframes: tuple[pd.DataFrame, pd.DataFrame]) -> HDXMeasurement:
    apo_corrected, _ = dataframes
    hdxm_apo = HDXMeasurement(apo_corrected, sequence=sequence, temperature=TEMPERATURE, pH=PH)

    return hdxm_apo


@pytest.fixture()
def hdxm_apo_red(dataframes: tuple[pd.DataFrame, pd.DataFrame]) -> HDXMeasurement:
    apo_corrected, _ = dataframes

    reduced_data = apo_corrected.copy()[apo_corrected["end"] < 40]

    hdxm_apo_red = HDXMeasurement(reduced_data, temperature=TEMPERATURE, pH=PH, c_term=155)

    return hdxm_apo_red


@pytest.fixture()
def hdxm_dimer(dataframes: tuple[pd.DataFrame, pd.DataFrame]) -> HDXMeasurement:
    _, dimer_corrected = dataframes
    hdxm_dimer = HDXMeasurement(
        dimer_corrected, sequence=sequence_dimer, temperature=TEMPERATURE, pH=PH
    )

    return hdxm_dimer


@pytest.fixture()
def hdxm_set(hdxm_apo: HDXMeasurement, hdxm_dimer: HDXMeasurement) -> HDXMeasurementSet:
    hdxm_set = HDXMeasurementSet([hdxm_apo, hdxm_dimer])
    return hdxm_set


# @pytest.fixture()
# def hdxm_apo_red(dataset: HDXDataSet) -> HDXMeasurement:
#     hdxm = HDXMeasurement.from_dataset(
#         dataset, state="SecB_tetramer_red", d_percentage=100.0, drop_first=1
#     )

#     return hdxm


# @pytest.fixture()
# def hdxm_set() -> HDXMeasurementSet:
#     yaml_file = input_dir / "data_states_deltas.yaml"
#     hdx_spec = yaml.safe_load(yaml_file.read_text())

#     dataset = HDXDataSet.from_spec(hdx_spec, data_dir=input_dir)
#     hdxm_set = HDXMeasurementSet.from_dataset(dataset)

#     return hdxm_set


@pytest.mark.skip(reason="Different result on py>3.9")
def test_initial_guess_wt_average(hdxm_apo_red: HDXMeasurement):
    result = fit_rates_weighted_average(hdxm_apo_red)
    output = result.output

    assert output.size == 100
    check_rates = csv_to_dataframe(output_dir / "ecSecB_reduced_guess.csv")
    pd.testing.assert_series_equal(check_rates["rate"], output["rate"])


def test_initial_guess_half_time_interpolate(hdxm_apo_red: HDXMeasurement):
    result = fit_rates_half_time_interpolate(hdxm_apo_red)
    assert isinstance(result, GenericFitResult)
    assert result.output.index.name == "r_number"
    assert result.output["rate"].mean() == pytest.approx(0.04343354509254464)

    # todo additional tests:
    #  result = fit_rates_half_time_interpolate()


@pytest.mark.skip(reason="Hangs on GitHub Actions on Ubuntu")
def test_dtype_cuda(hdxm_apo: HDXMeasurement):
    check_deltaG = csv_to_dataframe(output_dir / "ecSecB_torch_fit.csv")
    initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

    # cfg.set("fitting", "device", "cuda")

    with cfg.context({"fitting.device": "cuda"}):
        gibbs_guess = hdxm_apo.guess_deltaG(initial_rates["rate"]).to_numpy()

        if torch.cuda.is_available():
            fr_global = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=1000, r1=2)
            out_deltaG = fr_global.output
            for field in ["dG", "k_obs", "covariance"]:
                assert_series_equal(
                    check_deltaG[field],
                    out_deltaG[hdxm_apo.name, field],
                    rtol=0.01,
                    check_dtype=False,
                    check_names=False,
                )
        else:
            raise AssertionError("CUDA not available")
            # with pytest.raises(AssertionError, match=r".* CUDA .*"):
            #     fr_global = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=1000, r1=2)


def test_dtype_cpu(hdxm_apo: HDXMeasurement):
    initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")
    check_deltaG = csv_to_dataframe(output_dir / "ecSecB_torch_fit.csv")

    gibbs_guess = hdxm_apo.guess_deltaG(initial_rates["rate"])

    with cfg.context({"fitting.device": "cpu", "fitting.dtype": "float32"}):
        fr_global = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=1000, r1=2)
        dg = fr_global.model.dG
        assert dg.dtype == torch.float32

        out_deltaG = fr_global.output
        for field in ["dG", "k_obs"]:
            assert_series_equal(
                check_deltaG["SecB WT apo", field],
                out_deltaG[hdxm_apo.name, field],
                rtol=0.01,
                check_dtype=False,
                check_names=False,
            )


def test_duptake_fit(hdxm_apo: HDXMeasurement):
    fr = fit_d_uptake(hdxm_apo, r1=0.5, repeats=3, verbose=False)
    check_d_uptake = csv_to_dataframe(output_dir / "ecSecB_d_uptake.csv")

    np.allclose(check_d_uptake, fr.output)


def test_global_fit(hdxm_apo: HDXMeasurement):
    initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

    # t0 = time.time()  # Very crude benchmarks
    gibbs_guess = hdxm_apo.guess_deltaG(initial_rates["rate"])
    fr_global = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=1000, r1=2)
    # t1 = time.time()

    # assert t1 - t0 < 5  # Fails sometimes
    out_deltaG = fr_global.output
    check_deltaG = csv_to_dataframe(output_dir / "ecSecB_torch_fit.csv")

    for field in ["dG", "covariance", "k_obs"]:
        assert_series_equal(
            check_deltaG["SecB WT apo", field],
            out_deltaG[hdxm_apo.name, field],
            rtol=0.01,
            check_names=False,
        )

    errors = fr_global.get_squared_errors()
    assert errors.shape == (1, hdxm_apo.Np, hdxm_apo.Nt)


@pytest.mark.skip(reason="Longer fit is not checked by default due to long computation times")
def test_global_fit_extended(hdxm_apo: HDXMeasurement):
    check_deltaG = csv_to_dataframe(output_dir / "ecSecB_torch_fit_epochs_20000.csv")
    initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")
    gibbs_guess = hdxm_apo.guess_deltaG(initial_rates["rate"])

    t0 = time.time()  # Very crude benchmarks
    fr_global = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=20000, r1=2)
    t1 = time.time()

    assert t1 - t0 < 50
    out_deltaG = fr_global.output
    for field in ["dG", "k_obs", "covariance"]:
        assert_series_equal(
            check_deltaG["SecB WT apo", field],
            out_deltaG[hdxm_apo.name, field],
            rtol=0.01,
            check_dtype=False,
            check_names=False,
        )

    errors = fr_global.get_squared_errors()
    assert errors.shape == (1, hdxm_apo.Np, hdxm_apo.Nt)


@pytest.mark.skip(reason="Longer fit is not checked by default due to long computation times")
def test_global_fit_extended_cuda(hdxm_apo: HDXMeasurement):
    check_deltaG = csv_to_dataframe(output_dir / "ecSecB_torch_fit_epochs_20000.csv")
    initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")
    gibbs_guess = hdxm_apo.guess_deltaG(initial_rates["rate"])

    with cfg.context({"fitting.device": "cuda", "fitting.dtype": "float32"}):
        fr_global = fit_gibbs_global(hdxm_apo, gibbs_guess, epochs=20000, r1=2)
        out_deltaG = fr_global.output

        for field in ["dG", "k_obs"]:
            assert_series_equal(
                check_deltaG[hdxm_apo.name, field],
                out_deltaG[hdxm_apo.name, field],
                rtol=0.01,
                check_dtype=False,
            )


def test_batch_fit(hdxm_apo: HDXMeasurement, hdxm_dimer: HDXMeasurement, tmp_path):
    hdx_set = HDXMeasurementSet([hdxm_dimer, hdxm_apo])
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

    check_df = csv_to_dataframe(output_dir / "ecSecB_batch.csv")
    states = [("SecB_tetramer", "SecB WT apo"), ("SecB_dimer", "SecB his dimer apo")]

    for s_state, p_state in states:
        from pandas.testing import assert_series_equal

        result = output[p_state]["dG"]
        test = check_df[p_state]["dG"]

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

    gibbs_guess = hdx_set[0].guess_deltaG(guess["rate"])  # Guesses from first measurement
    aligned_result = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, r1=2, r2=5, epochs=1000)
    output = aligned_result.output
    check_df = csv_to_dataframe(output_dir / "ecSecB_batch_aligned.csv")
    states = ["SecB his dimer apo", "SecB WT apo"]

    for state in states:
        from pandas.testing import assert_series_equal

        result = output[state]["dG"]
        test = check_df[state]["dG"]

        assert_series_equal(result, test, rtol=0.1)


# batch fit on delta N/C tail dataset
@pytest.mark.skip(reason="Needs adapting of `data_states_delta.yaml` to run")
def test_batch_fit_delta(hdxm_set, tmp_path):
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
