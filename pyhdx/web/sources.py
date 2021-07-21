import param
from bokeh.models import ColumnDataSource
import numpy as np
from pyhdx.models import Protein
import pandas as pd
from lumen.util import get_dataframe_schema

#todo refactor module to models?

from lumen.sources import Source, cached_schema


class DataFrameSource(Source):

    tables = param.Dict({}, doc="Dictionary of tables in this Source")

    updated = param.Event()

    dropna = param.Boolean(True, doc='Remove rows of all NaN when adding / selecting / removing  dataframes')
    # def __init__(self, **params):
    #     pass
        # super().__init__(**params)
        # if self.df.columns.nlevels == 1:
        #     self.tables = [self.name]
        #     self.multiindex = False
        # elif self.df.columns.nlevels == 2:
        #     self.multiindex = True
        #     self.tables = [] # todo populate tables for multiindex
        # else:
        #     raise ValueError("Currently column multiindex beyond two levels is not supported")

    def remove_df(self, table, name, level):
        raise NotImplementedError('Removing datafarmes not implemented')

        self.updated = True

    def add_df(self, df, table, names=None):
        """
        #Todo method for adding a table to multindex source

        Parameters
        ----------
        df

        Returns
        -------

        """

        # todo check if df already present, update?
        # Todo check for name collisions between level and column names


        target_df = self.tables[table]
        df = df.copy()
        if target_df.columns.nlevels != df.columns.nlevels:
            if isinstance(names, str):
                names = [names]
            if len(names) != target_df.columns.nlevels - df.columns.nlevels:
                raise ValueError(f"Insufficient names provided to match target dataframe multindex level {df.columns.nlevels}")

            if df.columns.nlevels == 1:
                cols = ((column,) for column in df.columns)
            else:
                cols = df.columns
            tuples = tuple((*names, *tup) for tup in cols)
            new_index = pd.MultiIndex.from_tuples(tuples, names=target_df.columns.names)
            df.columns = new_index

        new_df = pd.concat([target_df, df], axis=1)
        if self.dropna:
            new_df = new_df.dropna(how='all')

        self.tables[table] = new_df
        #todo check for row indices

        self.updated = True

    def get(self, table, **query):
        df = self.tables[table]

        # This means querying a field with the same name as higher-levels columns is not possible
        # Todo check for name collisions between level and column names
        while df.columns.nlevels > 1:
            selected_col = query.pop(df.columns.names[0], False)
            if selected_col:
                df = df[selected_col]
                if self.dropna:
                    df = df.dropna(how='all')  # These subsets will have padded NaN rows. Remove?
            else:
                break

        dask = query.pop('__dask', False)
        df = self._filter_dataframe(df, **query)

        return df

    @cached_schema
    def get_schema(self, table=None):
        schemas = {}
        for name in self.tables:
            if table is not None and name != table:
                continue
            df = self.get(name)
            schemas[name] = get_dataframe_schema(df)['items']['properties']
        return schemas if table is None else schemas[table]

    def get_unique(self, table=None, field=None, **query):
        """Get unique values for specified tables and fields"""
        print('deprecation candidate')
        unique_values = {}
        for name in self.tables:
            if table is not None and name != table:
                continue
            df = self.get(name, **query)
            if field is not None:
                unique_values[name] = df[field].unique()
            else:
                unique_values[name] = {field_name: df[field_name].unique() for field_name in df.columns}

            return unique_values if table is None else unique_values[table]


class DataSource(param.Parameterized):
    tags = param.List(doc='List of tags to specify the type of data in the dataobject.')
    source = param.ClassSelector(ColumnDataSource, doc='ColumnDataSource object which is used for graphical display')
    renderer = param.String(default='line')
    default_color = param.Color(default='#0611d4')  # todo get default color from css?

    def __init__(self, input_data, **params):
        #update to lumen / pandas dataframes
        self.render_kwargs = {k: params.pop(k) for k in list(params.keys()) if k not in self.param}
        #todo currently this override colors in dic
        super(DataSource, self).__init__(**params)
        self.df = self.format_data(input_data)
        default_color = 'color' if 'color' in self.df.columns else self.default_color  #df columns does not work with multiindex dfs
        self.render_kwargs['color'] = self.render_kwargs.get('color', default_color)

        self.source = ColumnDataSource(self.df, name=self.name)

    def __getitem__(self, item):
        return self.source.data.__getitem__(item)

    @property
    def export_df(self):
        df = pd.DataFrame({k: v for k, v in self.source.data.items() if not k.startswith('__')})
        return df

    def format_data(self, input_data):
        #todo allow dataframes
        if isinstance(input_data, np.ndarray):
            return pd.DataFrame(input_data)
        elif isinstance(input_data, dict):
            return pd.DataFrame(input_data)
        elif isinstance(input_data, Protein):
            return input_data.df
        else:
            raise TypeError("Invalid input data type")

    def to_numpy(self):
        raise NotImplementedError('Converting to numpy rec array not implemented')

    @property
    def scalar_fields(self):
        """Returns a list of names of fields with scalar dtype"""
        return [name for name, data in self.source.data.items() if np.issubdtype(data.dtype, np.number)]

    @property
    def y(self):
        """:class:`~numpy.ndarray`: Array of y values"""
        if 'y' in self.render_kwargs:
            try:
                return self.source.data[self.render_kwargs['y']]
            except TypeError:
                return None #  'y' might be a list of y values
        elif 'y' in self.source.data:
            return self.source.data['y']
        else:
            return None

    @property
    def x(self):
        """:class:`~numpy.ndarray`: Array of x values"""
        if 'x' in self.render_kwargs:
            return self.source.data[self.render_kwargs['x']]
        elif 'x' in self.source.data:
            return self.source.data['x']
        else:
            return None

    def update(self, data_source_obj):
        """
        Update the data and source object

        Parameters
        ----------
        data_source

        Returns
        -------


        """

        self.source.data.update(**data_source_obj.source.data)

    def resolve_tags(self, tags):
        # format ['tag1', ('tag2a', 'tag2b') ] = tag1 OR (tag2a AND tag2b)

        for tag in tags:
            if isinstance(tag, str):
                bool = tag in self.tags
                # if tag in self.tags:
                #     return True
            else:
                bool = all(sub_tag in self.tags for sub_tag in tag)
            if bool:
                return True
        return False

    #def _resolve_tag(self, tag):


class MultiIndexDataSource(DataSource):
    pass