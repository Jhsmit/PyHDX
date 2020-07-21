import pytest
import os
from pyhdx import PeptideMeasurements, PeptideMasterTable, KineticsSeries
from pyhdx.fileIO import read_dynamx
from pyhdx.support import np_from_txt
from pyhdx.expfact.kint import calculate_kint_for_sequence, calculate_kint_per_residue
import numpy as np
from functools import reduce
from operator import add

directory = os.path.dirname(__file__)


class TestUptakeFileModels(object):

    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        cls.pf1 = PeptideMasterTable(read_dynamx(fpath))

        fpath = os.path.join(directory, 'test_data', 'ds2.csv')
        cls.pf2 = PeptideMasterTable(read_dynamx(fpath))

        fpath = os.path.join(directory, 'test_data', 'ds3.csv')
        cls.pf3 = PeptideMasterTable(read_dynamx(fpath))

    def test_peptidemastertable(self):
        data = self.pf1.data[self.pf1.data['start'] < 50]
        res = self.pf1.isin_by_idx(self.pf1.data, data)
        assert res.shape == self.pf1.data.shape
        assert len(data) == np.sum(res)

    def test_peptidemeasurement(self):
        assert isinstance(self.pf1, PeptideMasterTable)

        p_dict = self.pf1.return_by_name('PpiA-FD', 0.167)

        name = 'PpiANative_0.167'
        pm = p_dict[name]
        assert pm.name == name
        assert isinstance(pm, PeptideMeasurements)

        res_scores = np.arange(pm.prot_len)
        scores = pm.calc_scores(res_scores)
        assert len(scores) == len(pm)


    # def test_apply_controls(self):
    #     states_dict = self.pf1.groupby_state()
    #
    #     series1 = states_dict['PpiANative']
    #     assert len(series1) == 7
    #     assert len(series1.times) == 7
    #     assert isinstance(series1, KineticsSeries)
    #     control_100 = self.pf1.get_data('PpiA-FD', 0.167)
    #
    #     for pm in series1:
    #         assert len(pm) == 78
    #     series1.set_control(control_100)
    #     for pm in series1:
    #         assert len(pm) == 78

    def test_split(self):
        #control_100 = ("Full Deuteration control", 0.167)
        series_name = 'SecA-monomer'

        states = self.pf3.groupby_state()

        series = states[series_name]
        series.make_uniform()

        split_series = series.split()
        new_len = reduce(add, [reduce(add, [len(pm.data) for pm in ks]) for ks in split_series.values()])

        assert len(series.full_data) == new_len

        k1 = list(split_series.keys())[0]
        for k, v in split_series.items():
            s, e = np.array(k.split('_')).astype(int)
            pm = v[0]
            assert np.min(pm.data['start']) == s
            assert np.all(pm.data['end'] < e)

    def test_uniform(self):
        control_100 = ("Full Deuteration control", 0.167)
        series_name = 'SecA-monomer'

        states = self.pf3.groupby_state()
        series = states[series_name]

        assert ~series.uniform
        # new_series = series.make_uniform(in_place=False)
        # assert new_series.uniform
        series.make_uniform()
        assert series.uniform

        l1 = len(series[0])
        series.make_uniform()
        l2 = len(series[0])
        assert l1 == l2

    def test_ds2(self):
        self.pf2.set_control(('FD', 0.001), ('Native folded', 60.000004), remove_nan=False)
        states = self.pf2.groupby_state()
        series = states['folding_4C_10secLabelling']
        assert len(series[0]) == 80
        assert np.all(np.isnan(series[0].scores_average))

        self.pf2.set_control(('FD', 0.001), ('Native folded', 60.000004), remove_nan=True)
        states = self.pf2.groupby_state()
        series = states['folding_4C_10secLabelling']
        assert ~np.all(np.isnan(series[0].scores_average))
        assert len(series[0]) == 79


class TestSeries(object):
    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'ds1.csv')
        cls.pf1 = PeptideMasterTable(read_dynamx(fpath))

    def test_coverage(self):
        states = self.pf1.groupby_state()
        series = states['PpiANative']
        sequence = 'MFKSTLAAMAAVFALSALSPAAMAAKGDPHVLLTTSAGNIELELDKQKAPVSVQNFVDYVNSGFYNNTTFHRVIPGFMIQGGGFTEQMQQKKPNPPIKNEADNGLRNTRGTIAMARTADKDSATSQFFINVADNAFLDHGQRDFGYAVFGKVVKGMDVADKISQVPTHDVGPYQNVPSKPVVILSAKVLP'
        k_full = calculate_kint_for_sequence(1, len(sequence), sequence, 300, 7)
        k_part = series.cov.calc_kint(300, 7, None)

        for s1, s2, k1, k2 in zip(sequence, series.cov.sequence, k_full, k_part):
            if s2 == 'X':
                continue
            else:
                assert s1 == s2
                if s1 == 'P':
                    continue
                assert k1 == k2


class TestSimulatedData(object):
    @classmethod
    def setup_class(cls):
        fpath = os.path.join(directory, 'test_data', 'simulated_data.csv')
        cls.data = np_from_txt(fpath, delimiter=',')
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
        assert series.uniform
        assert len(series) == len(self.timepoints)
        peptides = series[3]

        assert peptides.start == self.start
        assert peptides.end == self.end
        assert len(peptides.r_number) == self.end - self.start
        assert np.all(np.diff(peptides.r_number) == 1)
        assert peptides.prot_len == self.end - self.start

        blocks = [1, 4, 2, 4, 3, 2, 10, 4, 2, 3, 3, 2, 1]
        assert np.all(blocks == peptides.block_length)

        lengths = peptides.data['end'] - peptides.data['start']
        assert np.all(lengths == peptides.data['ex_residues'])

        block_coverage = [True, True, True, True, True, True, True, False, True, True, True, True, True]
        assert np.all(block_coverage == peptides.block_coverage)

        for row in peptides.X:
            assert np.sum(row) == 1
        assert peptides.X.shape == (len(self.data) / len(self.timepoints), self.end - self.start)

        assert np.all(np.sum(peptides.X, axis=1) == 1)

        for row, elem in zip(peptides.X_red, peptides.data):
            assert np.nansum(row) == len(elem['sequence'])

        assert peptides.X_red.shape == (len(self.data) / len(self.timepoints), len(blocks))
        unique_elems = [np.unique(col[col != 0]) for col in peptides.X_red.T.astype(int)]
        cov_blocks = [elem[0] for elem in unique_elems if len(elem) == 1]
        assert np.all(peptides.block_length[peptides.block_coverage] == cov_blocks)

        assert peptides.block_length[~peptides.block_coverage] == self.nc_end - self.nc_start
        assert np.sum(peptides.has_coverage) == self.nc_start - self.start + self.end - self.nc_end

        assert peptides.exposure == self.timepoints[3]
        assert peptides.state == 'state1'
        assert peptides.sequence == self.sequence

        # series keys are inclusive, exclusive
        keys = [f'{self.start}_{self.nc_start}', f'{self.nc_end}_{self.end}']
        split = series.split()

        for k1, k2 in zip(keys, split.keys()):
            assert k1 == k2

        s1 = split[keys[0]]
        p1 = s1[3]
        assert p1.start == self.start
        assert p1.end == self.nc_start
        assert np.all(p1.r_number == np.arange(self.start, self.nc_start))

        s2 = split[keys[1]]
        p2 = s2[3]
        assert p2.start == self.nc_end
        assert p2.end == self.end
        assert np.all(p2.r_number == np.arange(self.nc_end, self.end))

        for i, t in enumerate(self.timepoints):
            # Check all timepoints in both split series
            assert s1[i].exposure == t
            assert s2[i].exposure == t

            total_peptides = np.sum([len(v[i].data) for v in split.values()])
            assert total_peptides == len(peptides.data)

    def test_drop_first_prolines(self):
        for i, df in enumerate([1, 2, 3]):
            print('df', df)
            pcf = PeptideMasterTable(self.data, drop_first=df, ignore_prolines=True, remove_nan=False)
            states = pcf.groupby_state()
            assert len(states) == 1
            series = states['state1']
            assert series.uniform
            assert len(series) == len(self.timepoints)

            peptides = series[3]
        #    assert peptides.start == self.start + df + 2 # 2 prolines
            assert peptides.end == self.end

            #This still includes peptides
            assert peptides.sequence[peptides.start:].lstrip('X') == self.sequence[peptides.start:]

            #take only the exchanging residues
            ex_res = ''.join(list(peptides.sequence[i] for i in peptides.r_number - 1))
            # this only holds up to the first coverage break
            assert ex_res[:10] == self.sequence[peptides.start - 1:].replace('P', '')[:10]


            assert np.sum(peptides.block_length) == len(peptides.r_number)
       #     assert len(peptides.r_number) == self.end - self.start + 1 - df  - 2 # 2 prolines

#            assert np.all(np.diff(peptides.r_number) == 1)  #THIS WILL CHANGE!


            #assert peptides.prot_len == self.end - self.start + 1 - df

            # Number of exchanging residues removed from peptides by drop_first and prolines
            reductions = [
                [4, 3, 2, 3, 2, 2, 1],
                [4, 3, 3, 4, 3, 2, 2],
                [4, 4, 4, 5, 4, 3, 3]
            ][i]

            print(i)
            #unmodified: [11 15  9 19  8  9  5]
            lengths = np.array([len(seq) for seq in peptides.data['sequence']]) - np.array(reductions)
            assert np.all(lengths == peptides.data['ex_residues'])

            # blocks = [
            #     [1,  4,  2,  3,  3,  2, 10,  5,  2,  3,  2,  2,  1],
            #     [1,  4,  2,  2,  3,  2, 10,  6,  2,  3,  1,  2,  1],
            #     [1,  4,  2,  1,  3,  2, 10,  7,  2,  3,  2,  1]
            # ][i]
            # assert np.all(blocks == peptides.block_length)
            #
            # # The block after the 10 length blocks should be the gap
            # assert not peptides.block_coverage[blocks.index(10) + 1]
            #
            # for row in peptides.X:
            #     assert np.round(np.sum(row), 6) == 1.
            # assert peptides.X.shape == (len(self.data) / len(self.timepoints), self.end - self.start + 1 - df)

