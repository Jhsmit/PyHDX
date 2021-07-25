from lumen.filters import Filter
from lumen.sources import Source
from param.parameterized import default_label_formatter
import param
from pyhdx.web.widgets import ColoredStaticText
import panel as pn


class WebAppFilter(Filter):
    """

    """

    source = param.ClassSelector(Source)

    # maybe instead of making filters co-dependent the second filter should have a DerivedSource
    # but we'll deal with this later
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

    @param.depends('value', watch=True)
    def update(self, *events):
        """gets called when the source event triggers"""
        self.updated = True


class WebAppWidgetFilter(WebAppFilter):

    empty_select = param.Boolean(default=True)  # currently unused param

    value = param.Selector()

    #_widget not defined: WebAppwidget is abstract base class
    # also upsdat is not defined

    def __init__(self, **params):
        super().__init__(**params)
        name = default_label_formatter(self.field)
        self.widget = self._widget.from_param(self.param.value, name=name)
        self.update()  # populate options

    @property
    def query(self):
        return self.widget.value

    @property
    def panel(self):
        return self.widget


class UniqueValuesFilter(WebAppWidgetFilter):
    """
    Selects a column from the data specified by 'fierld' and returns a select with options the unique values in this
    column

    """

    filter_type = 'unique_values'
    show_index = param.Boolean(False, doc='Set True to display the index of the unique values rather than their value')
    #_widget = pn.widgets.DiscreteSlider  #todo change back to DiscreteSlider, deal with initial value
    # see perhaps: https://github.com/holoviz/panel/pull/1837 ?
    _widget = pn.widgets.Select

    def __init__(self, **params):
        super().__init__(**params)

        self.unique_vals = []

    def update(self, *events):
        data = self.get_data()
        data = data.dropna(how='all') #todo perhaps add this line in get_data method of DataFrameSource

        if self.field not in data: # no data loaded
            options = []
        else:
            self.unique_vals = list(data[self.field].unique())

            if self.show_index:
                options = list(range(len(self.unique_vals)))
            else:
                options = self.unique_vals

        self.param['value'].objects = options
        if not self.value and options:
            self.value = options[0]

        self.updated = True

    @property
    def query(self):
        if self.show_index:
            try:  # try/except to handle initation with empty data
                return self.unique_vals[self.widget.value]
            except TypeError:
                return None
        else:
            return self.widget.value


class MultiIndexSelectFilter(WebAppWidgetFilter):
    """
    Select sub-dataframes from column-multiindex dataframes by their top-level index

    Can be chained together to select through multiple levels of dataframes

    filter1 = MultiIndexSelectFilter(source=source, table='my_table', field='name_top_level_indices')
    filter2 = MultiIndexSelectFilter(source=source, table='my_table', field='name_second_level_indices' filters=[filter1])
    filter3 = MultiIndexSelectFilter(source=source, table='my_table', field='name_third_level_indices', filters=[filter1, filter2])

    and SomeView(... filters=[filter1, filter2, filter3]

    #todo Make stacked filter into a single class


    """
    #todo subclass with uniquevaluesfilter
    # only create update param connection with source
    # needs some kind of method that returns the df its supposed to filter on (resolving stacked filters)

    filter_type = 'select column'

    # wildcard = param.Boolean(False, doc="Add a wildcard entry (*) to this filter's options")

    _widget = pn.widgets.Select

    #todo allow for none/wildcard (for losses)

    def __init__(self, **params):
        super().__init__(**params)
        #todo this is already done in superclass
        for filter in self.filters:
            filter.param.watch(self.update, 'updated')

    def update(self, *events):
        data = self.get_data()
        options = list(data.columns.get_level_values(0).unique())

        self.options = options  #does this update the widget?
        self.param['value'].objects = options
        if not self.value and options:
            self.value = options[0]

        self.updated = True

