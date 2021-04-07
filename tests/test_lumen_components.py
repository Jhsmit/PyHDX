from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable
import pandas as pd


from pyhdx.panel.data_sources import DataFrameSource
from pyhdx.panel.filters import MultiIndexSelectFilter

from lumen.sources import DerivedSource


directory = Path(__file__).parent

data_dir = directory / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))
states = pmt.groupby_state()

st1 = states['SecB his dimer apo']
st2 = states['SecB WT apo']

df1 = pd.DataFrame(st1.full_data)
df2 = pd.DataFrame(st2.full_data)


class TestLumenSources(object):
    @classmethod
    def setup_class(cls):
        cls.df1 = df1
        cls.df2 = df2


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

