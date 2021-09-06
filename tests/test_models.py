import pytest
import os
from pyhdx import PeptideMeasurements, PeptideMasterTable, HDXMeasurement
from pyhdx.models import Protein, Coverage
from pyhdx.fileIO import read_dynamx, csv_to_protein, csv_to_hdxm
import numpy as np
from functools import reduce
from operator import add
from pathlib import Path
import pandas as pd
import tempfile
import pickle
import pytest


directory = Path(__file__).parent


class TestUptakeFileModels(object):

    @classmethod
    def setup_class(cls):
        fpath = directory / 'test_data' / 'ecSecB_apo.csv'
        cls.pf1 = PeptideMasterTable(read_dynamx(fpath))

    def test_peptidemastertable(self):
        data = self.pf1.data[self.pf1.data['start'] < 50]
        res = self.pf1.isin_by_idx(self.pf1.data, data)
        assert res.shape == self.pf1.data.shape
        assert len(data) == np.sum(res)

    def test_peptidemeasurement(self):
        assert isinstance(self.pf1, PeptideMasterTable)

        self.pf1.set_control(('Full deuteration control', 0.167*60))

        data = self.pf1.get_state('SecB WT apo')
        assert isinstance(data, np.ndarray)


class TestHDXMeasurement(object):
    @classmethod
    def setup_class(cls):
        fpath = directory / 'test_data' / 'ecSecB_apo.csv'
        cls.pmt = PeptideMasterTable(read_dynamx(fpath))
        cls.pmt.set_control(('Full deuteration control', 0.167*60))
        d = cls.pmt.get_state('SecB WT apo')
        cls.temperature, cls.pH = 273.15 + 30, 8.
        cls.hdxm = HDXMeasurement(d, temperature=cls.temperature, pH=cls.pH)

    def test_dim(self):
        assert self.hdxm.Nt == len(np.unique(self.hdxm.full_data['exposure']))

    def test_guess(self):
        pass

    def test_tensors(self):
        tensors = self.hdxm.get_tensors()

        # assert ...

    def test_to_file(self):
        with tempfile.TemporaryDirectory() as tempdir:
            fpath = Path(tempdir) / 'hdxm.csv'
            self.hdxm.to_file(fpath)
            hdxm_read = csv_to_hdxm(fpath)
            k1 = self.hdxm.coverage['k_int']
            k2 = hdxm_read.coverage['k_int']
            pd.testing.assert_series_equal(k1, k2)

            assert self.hdxm.metadata == hdxm_read.metadata

@pytest.mark.skip(reason="Simulated data was removed")
class TestSimulatedData(object):
    @classmethod
    def setup_class(cls):
        fpath = directory / 'test_data' / 'simulated_data.csv'
        cls.data = txt_to_np(fpath, delimiter=',')
        cls.data['end'] += 1 # because this simulated data is in old format of inclusive, inclusive
        cls.sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

        cls.timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
        cls.start, cls.end = 5, 46  # total span of protein (inc, excl)
        cls.nc_start, cls.nc_end = 31, 35  # span of no coverage area (inc, inc)

    def test_keep_prolines(self):
        pcf = PeptideMasterTable(self.data, drop_first=0, ignore_prolines=False, remove_nan=False)
        states = pcf.groupby_state()
        assert len(states) == 1
        series = states['state1']
        assert len(series) == len(self.timepoints)
        peptides = series[3]

        assert len(peptides.r_number) == self.end - self.start
        assert np.all(np.diff(peptides.r_number) == 1)

        blocks = [1, 4, 2, 4, 3, 2, 10, 4, 2, 3, 3, 2, 1]
        assert np.all(blocks == peptides.block_length)

        lengths = peptides.data['end'] - peptides.data['start']
        assert np.all(lengths == peptides.data['ex_residues'])

        # for row in peptides.X:
        #     assert np.sum(row) == 1

        assert peptides.X.shape == (len(self.data) / len(self.timepoints), self.end - self.start)

        #assert np.all(np.sum(peptides.X, axis=1) == 1)

        for row, elem in zip(peptides.X, peptides.data):
            assert np.nansum(row) == len(elem['sequence'])

        assert np.sum(peptides.protein['coverage']) == self.nc_start - self.start + self.end - self.nc_end
        assert peptides.exposure == self.timepoints[3]
        assert peptides.state == 'state1'
        assert ''.join(peptides.protein['sequence']) == self.sequence

    def test_drop_first_prolines(self):
        for i, df in enumerate([1, 2, 3]):
            pcf = PeptideMasterTable(self.data, drop_first=df, ignore_prolines=True, remove_nan=False)
            states = pcf.groupby_state()
            assert len(states) == 1
            series = states['state1']
            assert len(series) == len(self.timepoints)

            peptides = series[3]
            #assert peptides.start == self.start + df + 2 # 2 prolines
            #assert peptides.end == self.end

            #take only the exchanging residues
            # this only holds up to the first coverage break

            assert np.sum(peptides.block_length) == len(peptides.r_number)
            reductions = [
                [4, 3, 2, 3, 2, 2, 1],
                [4, 3, 3, 4, 3, 2, 2],
                [4, 4, 4, 5, 4, 3, 3]
            ][i]

            #unmodified: [11 15  9 19  8  9  5]
            lengths = np.array([len(seq) for seq in peptides.data['sequence']]) - np.array(reductions)
            assert np.all(lengths == peptides.data['ex_residues'])


class TestCoverage(object):
    @classmethod
    def setup_class(cls):
        fpath = directory / 'test_data' / 'ecSecB_apo.csv'
        cls.pmt = PeptideMasterTable(read_dynamx(fpath))
        d = cls.pmt.groupby_state()
        cls.series = d['SecB WT apo']
        cls.sequence = 'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQM' \
                       'AHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA'

    def test_sequence(self):
        data = self.series[0].data
        cov = Coverage(data)

        for r, s in zip(cov.r_number, cov['sequence']):
            if s != 'X':
                assert self.sequence[r - 1] == s

        assert cov.protein.c_term == 155
        cov_seq = Coverage(data, sequence=self.sequence)
        assert cov_seq.protein.c_term == len(self.sequence)

        for r, s in zip(cov_seq.r_number, cov_seq['sequence']):
            assert self.sequence[r - 1] == s

    def test_dim(self):
        cov = self.series.coverage
        assert cov.Np == len(np.unique(cov.data['sequence']))
        assert cov.Nr == len(cov.r_number)


class TestProtein(object):
    @classmethod
    def setup_class(cls):
        dtype = [('r_number', int), ('apple', float)]
        array1 = np.empty(15, dtype=dtype)
        array1['r_number'] = np.arange(15) + 3
        array1['apple'] = np.ones(15) * 12
        cls.array1 = array1

        dtype = [('r_number', int), ('apple', float), ('grapes', float)]
        array2 = np.empty(17, dtype=dtype)
        array2['r_number'] = np.arange(17) + 6
        array2['apple'] = np.ones(17) * 10
        array2['grapes'] = np.ones(17) * 15 + np.random.rand(17)
        cls.array2 = array2

        dtype = [('r_number', int), ('pear', float), ('banana', float)]
        array3 = np.empty(10, dtype=dtype)
        array3['r_number'] = np.arange(10) + 1
        array3['pear'] = np.random.rand(10) + 20
        array3['banana'] = -(np.random.rand(10) + 20)
        cls.array3 = array3
        metadata = {'temperature': 273.15, 'pH': 7.5, 'mutations': ['V123S', 'P234S']}
        cls.protein = csv_to_protein(directory / 'test_data' / 'ecSecB_info.csv')
        cls.protein.metadata = metadata

        fpath = directory / 'test_data' / 'ecSecB_apo.csv'
        pf1 = PeptideMasterTable(read_dynamx(fpath))
        #states = pf1.groupby_state(c_term=200)
        cls.series = HDXMeasurement(pf1.get_state('SecB WT apo'), c_term=200)

    def test_artithmetic(self):
        p1 = Protein(self.array1, index='r_number')
        p2 = Protein(self.array2, index='r_number')

        comparison_quantity = 'apple'
        datasets = [self.array1, self.array2]  # len 15, 3 tm 17 , len 17, 6 tm 22
        r_all = np.concatenate([array['r_number'] for array in datasets])
        r_full = np.arange(r_all.min(), r_all.max() + 1)

        output = np.full_like(r_full, fill_value=np.nan,
                              dtype=[('r_number', int), ('value1', float), ('value2', float), ('comparison', float)])

        output['r_number'] = r_full
        idx = np.searchsorted(output['r_number'], self.array1['r_number'])
        output['value1'][idx] = self.array1[comparison_quantity]

        idx = np.searchsorted(output['r_number'], self.array2['r_number'])
        output['value2'][idx] = self.array2[comparison_quantity]

        comparison = p1 - p2
        output['comparison'] = output['value1'] - output['value2']
        assert np.allclose(comparison['apple'], output['comparison'], equal_nan=True)

        comparison = p1 + p2
        output['comparison'] = output['value1'] + output['value2']
        assert np.allclose(comparison['apple'], output['comparison'], equal_nan=True)

        comparison = p1 / p2
        output['comparison'] = output['value1'] / output['value2']
        assert np.allclose(comparison['apple'], output['comparison'], equal_nan=True)

        comparison = p1 * p2
        output['comparison'] = output['value1'] * output['value2']
        assert np.allclose(comparison['apple'], output['comparison'], equal_nan=True)

        data = np.random.rand(len(p1))
        pd_series = pd.Series(data, index=p1.index, name='apple')
        sub_result = p1.sub(pd_series, axis='index')
        assert np.allclose(sub_result['apple'], p1['apple'] - data)
        assert isinstance(sub_result, Protein)

    def test_k_int(self):
        protein = self.protein.copy()

        protein.df.rename(columns={'k_int': 'k_int_saved'}, inplace=True)
        protein.set_k_int(273.15 + 30, 8.)

        assert np.allclose(protein['k_int'], protein['k_int_saved'])

        # ecSecB
        self.series.coverage.protein.set_k_int(300., 8.)

        k_int = self.series.coverage.protein['k_int'].to_numpy()
        assert k_int[0] == np.inf  # N terminal exchange rate is zero
        assert np.all(k_int[-10:] == 0.)
        assert len(k_int) == self.series.coverage.protein.c_term

        prolines = self.series.coverage.protein['sequence'].to_numpy() == 'P'
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
            fpath = Path(tempdir) / 'protein.csv'
            self.protein.to_file(fpath)
            protein_read = csv_to_protein(fpath)
            pd.testing.assert_frame_equal(self.protein.df, protein_read.df)
            assert self.protein.metadata == protein_read.metadata
