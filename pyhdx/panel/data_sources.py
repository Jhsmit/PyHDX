import param
from bokeh.models import ColumnDataSource


class DataSource(param.Parameterized):
    tags = param.List()
    source = param.ClassSelector(ColumnDataSource)
    renderer = param.String(default='line')

    def __init__(self, dic, **params):
        super(DataSource, self).__init__(**params)
        self.source = ColumnDataSource(dic)




ds = DataSource({}, tags=['123', 'asdf'])
print(ds.tags)