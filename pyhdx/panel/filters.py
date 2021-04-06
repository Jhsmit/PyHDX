from lumen.filters import Filter
from lumen.sources import Source
from param.parameterized import default_label_formatter
import param
from pyhdx.panel.widgets import ColoredStaticText
import panel as pn


class WebAppFilter(Filter):
    """
    
    """

    source = param.ClassSelector(Source)

    filters = param.List()
    updated = param.Event()

    def __init__(self, **params):
        super().__init__(**params)

        if self.source:
            self.source.param.watch(self.update, 'updated')

        for filt in self.filters:
            filt.param.watch(self.update, 'updated')

    def get_data(self):
        """Equivalent to view's get_data method"""

        query = {
            filt.field: filt.query for filt in self.filters
            if filt.query is not None and
            (filt.table is None or filt.table == self.table)
        }

        data = self.source.get(self.table, **query)
        return data

    def update(self, *events):
        """gets called when the source event triggers"""
        self.updated = True


class WebAppWidgetFilter(WebAppFilter):

    #options = param.List(default=[], doc="""Unique values in specified field to filter""")

    empty_select = param.Boolean(default=True)

    value = param.Selector()

    #_widget not defined: WebAppwidget is abstract base class
    # also upsdat is not defined

    def __init__(self, **params):
        super().__init__(**params)
        name = default_label_formatter(self.field)
        self.widget = self._widget.from_param(self.param.value, name=name)
        self.update()  # populate options

    @param.depends('value', watch=True)
    def _print(self):
        # is there another way of linking these?
        # parameter js_link?
        self.updated = True

    @property
    def query(self):
        return self.widget.value

    @property
    def panel(self):
        return self.widget


class UniqueValuesFilter(WebAppWidgetFilter):

    filter_type = 'unique_values'
    #_widget = pn.widgets.DiscreteSlider  #todo change back to DiscreteSlider, deal with initial value
    _widget = pn.widgets.Select

    def __init__(self, **params):
        super().__init__(**params)

    def update(self, *events):
        data = self.get_data()
        data = data.dropna(how='all')  #todo perhaps add this line in get_data method of DataFrameSource
        options = list(data[self.field].unique())

        self.param['value'].objects = options
        if not self.value:
            self.value = options[0]

        self.updated = True


class SelectFilter(WebAppWidgetFilter):
    #todo subclass with uniquevaluesfilter
    # only create update param connection with source
    # needs some kind of method that returns the df its supposed to filter on (resolving stacked filters)

    filter_type = 'select column'

    _widget = pn.widgets.Select

    def __init__(self, **params):
        super().__init__(**params)

    def update(self, *events):
        data = self.get_data()
        options = list(data.columns.levels[0])

        self.options = options  #does this update the widget?
        self.param['value'].objects = options
        if not self.value:
            self.value = options[0]

        self.updated = True

