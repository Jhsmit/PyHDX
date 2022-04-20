import pyhdx
from pyhdx.fileIO import (
    read_dynamx,
    csv_to_dataframe,
    csv_to_protein,
    dataframe_to_stringio,
    dataframe_to_file,
    save_fitresult,
    load_fitresult,
)
from pyhdx.models import Protein, PeptideMasterTable, HDXMeasurement, HDXMeasurementSet
from pyhdx.fitting import fit_gibbs_global
from pathlib import Path
from io import StringIO
import numpy as np
import pandas as pd
import pytest

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"


class TestFileIO(object):
    @classmethod
    def setup_class(cls):
        cls.fpath = input_dir / "ecSecB_apo.csv"
        data = read_dynamx(cls.fpath)
        control = ("Full deuteration control", 0.167 * 60)

        cls.temperature, cls.pH = 273.15 + 30, 8.0

        pf = PeptideMasterTable(
            data, drop_first=1, ignore_prolines=True, remove_nan=False
        )
        pf.set_control(control)
        cls.hdxm = HDXMeasurement(
            pf.get_state("SecB WT apo"), temperature=cls.temperature, pH=cls.pH
        )

        initial_rates = csv_to_dataframe(output_dir / "ecSecB_guess.csv")

        gibbs_guess = cls.hdxm.guess_deltaG(initial_rates["rate"])
        cls.fit_result = fit_gibbs_global(cls.hdxm, gibbs_guess, epochs=100, r1=2)

    def test_read_dynamx(self):
        df = read_dynamx(self.fpath)

        assert df.shape[0] == 567
        assert df["start"][0] == 9
        assert df["end"][0] == 18
        df = read_dynamx(self.fpath, intervals=("exclusive", "inclusive"))
        assert df["start"][0] == 10
        assert np.sum(df["exposure"]) == 441632.55024
        assert np.sum(df["uptake"]) == 1910.589614

        df = read_dynamx(self.fpath, intervals=("inclusive", "exclusive"))
        assert df["end"][0] == 17

        with pytest.raises(ValueError):
            df = read_dynamx(self.fpath, intervals=("foo", "inclusive"))
        with pytest.raises(ValueError):
            df = read_dynamx(self.fpath, intervals=("inclusive", "bar"))

        with open(self.fpath, mode="r") as f:
            df = read_dynamx(StringIO(f.read()))
            assert df.shape[0] == 567

    def test_read_write_tables(self, tmp_path):
        # Single-index columns
        df = pd.DataFrame(np.random.randn(25, 4), columns=list("ABCD"))
        df.index.name = "singlecolumnindex"

        sio = StringIO()
        dataframe_to_stringio(df, sio)
        sio.seek(0)
        df_read = csv_to_dataframe(sio)
        pd.testing.assert_frame_equal(df, df_read)

        fpath = Path(tmp_path) / "single_index.csv"
        dataframe_to_file(fpath, df)
        csv_to_dataframe(fpath)
        pd.testing.assert_frame_equal(df, df_read)

        # multi-index column
        cols = pd.MultiIndex.from_product([("a", "b"), ("x", "y")])
        df = pd.DataFrame(np.random.randn(25, 4), columns=cols)
        df.index.name = "multicolumnindex"

        sio = StringIO()
        dataframe_to_stringio(df, sio)
        sio.seek(0)
        df_read = csv_to_dataframe(sio)
        pd.testing.assert_frame_equal(df, df_read)

        fpath = Path(tmp_path) / "multi_index.csv"
        dataframe_to_file(fpath, df)
        df_read = csv_to_dataframe(fpath)
        pd.testing.assert_frame_equal(df, df_read)

        protein = csv_to_protein(fpath)
        assert protein.index.name == "r_number"
        assert isinstance(protein, Protein)

        metadata = {
            "instrumuent": "LCMS",
            "settings": {"pressure": "5 kPa", "temperature": "400K"},
        }

        df.attrs["metadata"] = metadata

        fpath = Path(tmp_path) / "multi_index_with_metadata.csv"
        dataframe_to_file(fpath, df)
        df_read = csv_to_dataframe(fpath)
        pd.testing.assert_frame_equal(df, df_read)

        assert df_read.attrs["metadata"] == metadata

        fpath = Path(tmp_path) / "multi_index_with_metadata.txt"
        dataframe_to_file(fpath, df, fmt="pprint", include_version=True)
        lines = Path(fpath).read_text().split("\n")
        assert len(lines) == 38
        assert lines[0].strip() == pyhdx.VERSION_STRING

    def test_load_save_fitresult(self, tmp_path):
        # todo missing read batch result test

        fpath = Path(tmp_path) / "fit_result_single.csv"
        self.fit_result.to_file(fpath)
        df = csv_to_dataframe(fpath)
        assert df.attrs["metadata"] == self.fit_result.metadata
        fit_result_dir = Path(tmp_path) / "fit_result"

        save_fitresult(fit_result_dir, self.fit_result, log_lines=["test123"])

        log_lines = Path(fit_result_dir / "log.txt").read_text().split("\n")
        assert log_lines[-1] == "test123"

        fit_result_loaded = load_fitresult(fit_result_dir)
        assert isinstance(fit_result_loaded.losses, pd.DataFrame)
        assert isinstance(fit_result_loaded.hdxm_set, HDXMeasurementSet)

        timepoints = np.linspace(0, 30 * 60, num=100)
        d_calc = fit_result_loaded(timepoints)
        assert d_calc.shape == (1, self.hdxm.Np, len(timepoints))

        losses = csv_to_dataframe(fit_result_dir / "losses.csv")
        fr_load_with_hdxm_and_losses = load_fitresult(fit_result_dir)
        assert len(fr_load_with_hdxm_and_losses.losses) == 100

        assert (
            fr_load_with_hdxm_and_losses.metadata["total_loss"] == losses.iloc[-1].sum()
        )
