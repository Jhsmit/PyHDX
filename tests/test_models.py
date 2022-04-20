import pytest
import os
from pyhdx import HDXTimepoint, PeptideMasterTable, HDXMeasurement
from pyhdx.models import Protein, Coverage
from pyhdx.fileIO import read_dynamx, csv_to_protein, csv_to_hdxm, csv_to_dataframe
import numpy as np
from functools import reduce
from operator import add
from pathlib import Path
import pandas as pd
from pandas.testing import assert_frame_equal
import tempfile
import pickle
import pytest


cwd = Path(__file__).parent
input_dir = cwd / "test_data" / "input"
output_dir = cwd / "test_data" / "output"


class TestPeptideMasterTable(object):
    @classmethod
    def setup_class(cls):
        fpath = input_dir / "ecSecB_apo.csv"
        cls.pmt = PeptideMasterTable(read_dynamx(fpath))

    def test_properties(self):
        states = ["Full deuteration control", "SecB WT apo"]
        assert list(self.pmt.states) == states

        exposures = np.array(
            [0.0, 10.020000000000001, 30.0, 60.0, 300.0, 600.0, 6000.00048]
        )
        assert np.allclose(exposures, self.pmt.exposures, atol=0.01)

    # def test_methods(self):
    # testing of PeptideMasterTable methods


class TestHDXMeasurement(object):
    @classmethod
    def setup_class(cls):
        fpath = input_dir / "ecSecB_apo.csv"
        cls.pmt = PeptideMasterTable(read_dynamx(fpath))
        cls.pmt.set_control(("Full deuteration control", 0.167 * 60))
        d = cls.pmt.get_state("SecB WT apo")
        cls.temperature, cls.pH = 273.15 + 30, 8.0
        cls.hdxm = HDXMeasurement(d, temperature=cls.temperature, pH=cls.pH)

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
        cls.pmt = PeptideMasterTable(read_dynamx(fpath))
        data = cls.pmt.get_state("SecB WT apo")
        cls.hdxm = HDXMeasurement(data)
        cls.sequence = (
            "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQM"
            "AHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"
        )

    def test_sequence(self):
        data = self.hdxm[0].data
        cov = Coverage(data)

        for r, s in zip(cov.r_number, cov["sequence"]):
            if s != "X":
                assert self.sequence[r - 1] == s

        assert cov.protein.c_term == 155
        cov_seq = Coverage(data, sequence=self.sequence)
        assert cov_seq.protein.c_term == len(self.sequence)

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


class TestProtein(object):
    @classmethod
    def setup_class(cls):
        dtype = [("r_number", int), ("apple", float)]
        array1 = np.empty(15, dtype=dtype)
        array1["r_number"] = np.arange(15) + 3
        array1["apple"] = np.ones(15) * 12
        cls.array1 = array1

        dtype = [("r_number", int), ("apple", float), ("grapes", float)]
        array2 = np.empty(17, dtype=dtype)
        array2["r_number"] = np.arange(17) + 6
        array2["apple"] = np.ones(17) * 10
        array2["grapes"] = np.ones(17) * 15 + np.random.rand(17)
        cls.array2 = array2

        dtype = [("r_number", int), ("pear", float), ("banana", float)]
        array3 = np.empty(10, dtype=dtype)
        array3["r_number"] = np.arange(10) + 1
        array3["pear"] = np.random.rand(10) + 20
        array3["banana"] = -(np.random.rand(10) + 20)
        cls.array3 = array3
        metadata = {"temperature": 273.15, "pH": 7.5, "mutations": ["V123S", "P234S"]}
        cls.protein = csv_to_protein(output_dir / "ecSecB_info.csv")
        cls.protein.metadata = metadata

        fpath = input_dir / "ecSecB_apo.csv"
        pf1 = PeptideMasterTable(read_dynamx(fpath))
        cls.series = HDXMeasurement(pf1.get_state("SecB WT apo"), c_term=200)

    def test_artithmetic(self):
        p1 = Protein(self.array1, index="r_number")
        p2 = Protein(self.array2, index="r_number")

        comparison_quantity = "apple"
        datasets = [self.array1, self.array2]  # len 15, 3 tm 17 , len 17, 6 tm 22
        r_all = np.concatenate([array["r_number"] for array in datasets])
        r_full = np.arange(r_all.min(), r_all.max() + 1)

        output = np.full_like(
            r_full,
            fill_value=np.nan,
            dtype=[
                ("r_number", int),
                ("value1", float),
                ("value2", float),
                ("comparison", float),
            ],
        )

        output["r_number"] = r_full
        idx = np.searchsorted(output["r_number"], self.array1["r_number"])
        output["value1"][idx] = self.array1[comparison_quantity]

        idx = np.searchsorted(output["r_number"], self.array2["r_number"])
        output["value2"][idx] = self.array2[comparison_quantity]

        comparison = p1 - p2
        output["comparison"] = output["value1"] - output["value2"]
        assert np.allclose(comparison["apple"], output["comparison"], equal_nan=True)

        comparison = p1 + p2
        output["comparison"] = output["value1"] + output["value2"]
        assert np.allclose(comparison["apple"], output["comparison"], equal_nan=True)

        comparison = p1 / p2
        output["comparison"] = output["value1"] / output["value2"]
        assert np.allclose(comparison["apple"], output["comparison"], equal_nan=True)

        comparison = p1 * p2
        output["comparison"] = output["value1"] * output["value2"]
        assert np.allclose(comparison["apple"], output["comparison"], equal_nan=True)

        data = np.random.rand(len(p1))
        pd_series = pd.Series(data, index=p1.index, name="apple")
        sub_result = p1.sub(pd_series, axis="index")
        assert np.allclose(sub_result["apple"], p1["apple"] - data)
        assert isinstance(sub_result, Protein)

    def test_k_int(self):
        protein = self.protein.copy()

        protein.df.rename(columns={"k_int": "k_int_saved"}, inplace=True)
        k_int = protein.get_k_int(273.15 + 30, 8.0)

        assert np.allclose(k_int.to_numpy(), protein["k_int_saved"])

        # ecSecB
        k_int = self.series.coverage.protein.get_k_int(300.0, 8.0).to_numpy()

        assert k_int[0] == np.inf  # N terminal exchange rate is zero
        assert np.all(k_int[-10:] == 0.0)
        assert len(k_int) == self.series.coverage.protein.c_term

        prolines = self.series.coverage.protein["sequence"].to_numpy() == "P"
        assert np.all(k_int[prolines] == 0)

    def test_pickling(self):
        with tempfile.TemporaryFile() as tf:
            pickle.dump(self.protein, tf)
            tf.seek(0)
            unpickled = pickle.load(tf)

        pd.testing.assert_frame_equal(self.protein.df, unpickled.df)
        assert self.protein.metadata == unpickled.metadata

    def test_to_file(self):
        with tempfile.TemporaryDirectory() as tempdir:
            fpath = Path(tempdir) / "protein.csv"
            self.protein.to_file(fpath)
            protein_read = csv_to_protein(fpath)
            pd.testing.assert_frame_equal(self.protein.df, protein_read.df)
            assert self.protein.metadata == protein_read.metadata
