from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable
import pandas as pd
import numpy as np


from pyhdx.web.sources import DataFrameSource
from pyhdx.web.filters import MultiIndexSelectFilter

from lumen.sources import DerivedSource


directory = Path(__file__).parent

data_dir = directory / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167*60))

st1 = pmt.get_state('SecB his dimer apo')
st2 = pmt.get_state('SecB WT apo')

df1 = pd.DataFrame(st1)
df2 = pd.DataFrame(st2)

rates_df = pd.read_csv(data_dir / 'ecSecB_rates.txt', index_col=0, header=[0, 1])


class TestLumenSources(object):
    @classmethod
    def setup_class(cls):
        cls.df1 = df1
        cls.df2 = df2

    def test_adding_dataframes(self):
        #todo add test with non-alphabetical ordered columns (_deltaG, deltaG, covariance)
        col_index = pd.MultiIndex.from_tuples([], names=('fit ID', 'state', 'quantity'))
        row_index = pd.RangeIndex(0, 1, name='r_number')
        df_rates = pd.DataFrame(columns=col_index, index=row_index)

        tables = {'rates': df_rates}  # rates is nlevels == 3 dataframe
        source = DataFrameSource(tables=tables, name='dataframe', dropna=False)


        source.add_df(rates_df, 'rates', 'rates_fit')
        output_df = source.get('rates')

        assert output_df.size == 294
        level_1 = output_df.columns.levels[1]
        assert level_1.name == 'state'
        assert 'SecB WT apo' in level_1
        assert 'SecB his dimer apo' in level_1


        # Add nlvels == 1 dataframe
        df = pd.DataFrame({'rate': np.random.rand(100)})

        names = ['top_index', 'secb monomer']
        source.add_df(df, 'rates', names)
        output_df = source.get('rates')

        assert output_df.size == 468

        level_1 = output_df.columns.levels[1]
        assert level_1.name == 'state'
        assert 'secb monomer' in level_1

        level_0 = output_df.columns.levels[0]
        assert level_0.name == 'fit ID'
        assert 'top_index' in level_0


class TestLumenFilters(object):

    @classmethod
    def setup_class(cls):
        cls.df1 = df1
        cls.df2 = df2

        cls.df_level_2 = pd.concat([df1, df2], keys=['df1', 'df2'], names=['state', 'quantity'], axis=1)
        cls.df_level_3 = cls.df_level_2.copy()
        new_index = pd.MultiIndex.from_product([['top_index']] + cls.df_level_2.columns.levels,
                                               names=['top_name'] + cls.df_level_2.columns.names)

        cls.df_level_3.columns = new_index

        tables = {'level_1': cls.df1, 'level_2': cls.df_level_2, 'level_3': cls.df_level_3}
        cls.source = DataFrameSource(tables=tables)

    def test_multiindex_filters(self):
        filter = MultiIndexSelectFilter(source=self.source, table='level_2', field='state')
        assert filter.param['value'].objects == ['df1', 'df2']

        data = filter.get_data()
        assert data.columns.nlevels == 2

        ds = DerivedSource(source=self.source, filters=[filter])
        data = ds.get('level_2')
        assert data.columns.nlevels == 1

