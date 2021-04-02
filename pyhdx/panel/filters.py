from lumen.filters import Filter
from param.parameterized import default_label_formatter
import param
from pyhdx.panel.widgets import ColoredStaticText
import panel as pn

class UniqueValuesFilter(Filter):

    values = param.List(doc="""Unique values in specified field to filter""")

    filter_type = 'unique_values'

    def __init__(self, **params):
        super().__init__(**params)

        name = default_label_formatter(self.field)
        self.widget = pn.widgets.DiscreteSlider(name=name, options=self.values)
        self.widget.link(self, value='value')

    @property
    def query(self):
        return self.widget.value

    @property
    def panel(self):
        widget = self.widget.clone()
        self.widget.link(widget, value='value', bidirectional=True)
        return widget

    @param.depends('value', watch=True)
    def _print(self):
        print('filter value', self.value)

