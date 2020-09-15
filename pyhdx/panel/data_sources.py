import param
from bokeh.models import ColumnDataSource
import numpy as np


#todo refactor module to models?

class DataSource(param.Parameterized):
    tags = param.List(doc='List of tags to specify the type of data in the dataobject.')
    source = param.ClassSelector(ColumnDataSource, doc='ColumnDataSource object which is used for graphical display')
    renderer = param.String(default='line')
    color = param.Color(default='#0611d4')  # todo get default color from css?

    def __init__(self, input_data, **params):
        self.render_kwargs = {k: params.pop(k) for k in list(params.keys()) if k not in self.param}
        self.render_kwargs['color'] = self.render_kwargs.get('color', 'color')
        super(DataSource, self).__init__(**params)
        dic = self.get_dic(input_data)
        self.source = ColumnDataSource(dic)

    def get_dic(self, input_data):
        if isinstance(input_data, np.ndarray):
            dic = {name: input_data[name] for name in input_data.dtype.names}
            #self.array = input_data  #
        elif isinstance(input_data, dict):
            dic = input_data
            #self.array = None
        else:
            raise TypeError("Invalid input data type")

        if 'color' not in dic.keys():
            column = next(iter(dic.values()))
            color = np.full_like(column, fill_value=self.color, dtype='U7')
            dic['color'] = color
        return dic

    @property
    def y(self):
        """:class:`~numpy.ndarray`: Array of y values"""
        try:
            y_name = self.render_kwargs['y']
            return self.source.data[y_name]
        except KeyError:
            return None

    @property
    def x(self):
        """:class:`~numpy.ndarray`: Array of x values"""
        try:
            x_name = self.render_kwargs['x']
            return self.source.data[x_name]
        except KeyError:
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