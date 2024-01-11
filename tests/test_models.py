from pyhdx import HDXMeasurement
from pyhdx.datasets import read_dynamx
from pyhdx.models import Coverage
from pyhdx.fileIO import csv_to_hdxm, csv_to_dataframe
import numpy as np
from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal
import tempfile

from pyhdx.process import apply_control, correct_d_uptake, filter_peptides

cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"


class TestHDXMeasurement(object):
    @classmethod
    def setup_class(cls):
        fpath = input_dir / "ecSecB_apo.csv"
        df = read_dynamx(fpath)

        fd = {
            "state": "Full deuteration control",
            "exposure": {"value": 0.167, "unit": "min"},
        }

        fd_df = filter_peptides(df, **fd)
        peptides = filter_peptides(df, state="SecB WT apo")  # , query=["exposure != 0."])
        peptides_control = apply_control(peptides, fd_df)
        peptides_corrected = correct_d_uptake(peptides_control)

        cls.temperature, cls.pH = 273.15 + 30, 8.0
        cls.hdxm = HDXMeasurement(
            peptides_corrected, temperature=cls.temperature, pH=cls.pH, c_term=155
        )

    def test_dim(self):
        assert self.hdxm.Nt == len(self.hdxm.data["exposure"].unique())

    def test_guess(self):
        pass

    def test_tensors(self):
        tensors = self.hdxm.get_tensors()
        # assert ...

    def test_rfu(self):
        rfu_residues = self.hdxm.rfu_residues
        compare = csv_to_dataframe(output_dir / "ecSecB_rfu_per_exposure.csv")
        compare.columns = compare.columns.astype(float)
        compare.columns.name = "exposure"
        assert_frame_equal(rfu_residues, compare)

    def test_to_file(self):
        with tempfile.TemporaryDirectory() as tempdir:
            fpath = Path(tempdir) / "hdxm.csv"
            self.hdxm.to_file(fpath)
            hdxm_read = csv_to_hdxm(fpath)
            k1 = self.hdxm.coverage["k_int"]
            k2 = hdxm_read.coverage["k_int"]
            pd.testing.assert_series_equal(k1, k2)

            assert self.hdxm.metadata == hdxm_read.metadata


class TestCoverage(object):
    @classmethod
    def setup_class(cls):
        fpath = input_dir / "ecSecB_apo.csv"
        df = read_dynamx(fpath)

        fd = {
            "state": "Full deuteration control",
            "exposure": {"value": 0.167, "unit": "min"},
        }

        fd_df = filter_peptides(df, **fd)
        peptides = filter_peptides(df, state="SecB WT apo")  # , query=["exposure != 0."])
        peptides_control = apply_control(peptides, fd_df)
        peptides_corrected = correct_d_uptake(peptides_control)

        cls.hdxm = HDXMeasurement(peptides_corrected, c_term=155)
        cls.sequence = (
            "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQM"
            "AHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"
        )

    def test_sequence(self):
        data = self.hdxm[0].data
        cov = Coverage(data, c_term=155)

        for r, s in zip(cov.r_number, cov["sequence"]):
            if s != "X":
                assert self.sequence[r - 1] == s

        assert cov.protein.index.max() == 155
        cov_seq = Coverage(data, sequence=self.sequence)
        assert cov_seq.protein.index.max() == len(self.sequence)

        for r, s in zip(cov_seq.r_number, cov_seq["sequence"]):
            assert self.sequence[r - 1] == s

    def test_dim(self):
        cov = self.hdxm.coverage
        assert cov.Np == len(np.unique(cov.data["sequence"]))
        assert cov.Nr == len(cov.r_number)

        assert cov.Np == 63
        assert cov.Nr == 146

    def test_XZ(self):
        test_X = np.genfromtxt(output_dir / "attributes" / "X.txt")
        assert np.allclose(self.hdxm.coverage.X, test_X)

        test_Z = np.genfromtxt(output_dir / "attributes" / "Z.txt")
        assert np.allclose(self.hdxm.coverage.Z, test_Z)
