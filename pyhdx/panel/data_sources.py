import param
from bokeh.models import ColumnDataSource
import numpy as np
from pyhdx.models import Protein


#todo refactor module to models?

class DataSource(param.Parameterized):
    tags = param.List(doc='List of tags to specify the type of data in the dataobject.')
    source = param.ClassSelector(ColumnDataSource, doc='ColumnDataSource object which is used for graphical display')
    renderer = param.String(default='line')
    default_color = param.Color(default='#0611d4')  # todo get default color from css?

    def __init__(self, input_data, **params):
        self.render_kwargs = {k: params.pop(k) for k in list(params.keys()) if k not in self.param}
        #todo currently this override colors in dic
        super(DataSource, self).__init__(**params)
        dic = self.get_dic(input_data)
        default_color = 'color' if 'color' in dic else self.default_color
        self.render_kwargs['color'] = self.render_kwargs.get('color', default_color)

        self.source = ColumnDataSource(dic, name=self.name)

    def __getitem__(self, item):
        return self.source.data.__getitem__(item)

    def get_dic(self, input_data):
        if isinstance(input_data, np.ndarray):
            dic = {name: input_data[name] for name in input_data.dtype.names}
            #self.array = input_data  #
        elif isinstance(input_data, dict):
            dic = {k: np.array(v) for k, v in input_data.items()}
        elif isinstance(input_data, Protein):
            dic = {k: np.array(v) for k, v in input_data.to_dict('list').items()}
            dic['r_number'] = np.array(input_data.index)
        else:
            raise TypeError("Invalid input data type")

        #todo this does not apply to all data sets?
        if 'color' not in dic.keys():
            column = next(iter(dic.values()))
            color = np.full_like(column, fill_value=self.default_color, dtype='<U7')
            dic['color'] = color
        return dic

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