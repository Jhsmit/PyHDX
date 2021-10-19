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
    #  ^ this was as long time ago
    filters = param.List()

    transforms = param.List()

    updated = param.Event()

    pd_function = param.String(doc='Pandas function which this filter applies to the DataFrame ')

    def __init__(self, **params):
        super().__init__(**params)
        # self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        # super().__init__(**{k: v for k, v in params.items() if k in self.param})

        if self.source:
            self.source.param.watch(self.update, 'updated')

        for filt in self.filters:
            filt.param.watch(self.update, 'updated')

    def get_data(self):
        """Equivalent to view's get_data method"""

        queries = [filt.query for filt in self.filters]
        data = self.source.get(self.table, *queries)
        for transform in self.transforms:
            data = transform.apply(data)

        return data

    @param.depends('value', watch=True)
    def update(self, *events):
        """gets called when the source event triggers"""
        self.updated = True

    # def _get_params(self):
    #     return self.kwargs

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas functin to apply filter"""
        return {}

    @property
    def query(self):
        return self.pd_function, self.pd_kwargs

class TransformFilter(WebAppFilter):
    pd_function = param.String('transform')


class StackFilter(WebAppFilter):
    pd_function = 'stack'

    level = param.Integer(-1)  #actually int, str, or list of these, default -1 (last level)
    #level = param.ClassSelector(_class=[int, str, list])  # where list is list of (str | int)

    dropna = param.Boolean(True)

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        #todo get_params func which finds the correct params here
        return dict(level=self.level, dropna=self.dropna)

#todo some kind of class that can do multiple of these combined?

class YourMomFilter(WebAppFilter):
    pd_function = param.String()

    def __init__(self, **params):
        self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**{k: v for k, v in params.items() if k in self.param})

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas functin to apply filter"""
        return self.kwargs

class PivotFilter(WebAppFilter):
    pd_function = 'pivot'

    index = param.ClassSelector(class_=[str, list])

    columns = param.ClassSelector(class_=[str, list])

    values = param.ClassSelector(class_=[str, list])

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        #todo get_params func which finds the correct params here
        return dict(index=self.index, columns=self.columns, values=self.values)


class CrossSectionFilter(WebAppFilter):

    pd_function = 'xs'

    key = param.Tuple()

    axis = param.Integer(1, bounds=[0, 1])

    n_levels = param.Integer(None, doc="Number of levels")

    level = param.List()

    drop_level = param.Boolean(True)

    names = param.List(None, doc="List of label names for widgets")

    empty_select = param.Boolean(default=False, doc="""
        Add an option to Select widgets to indicate no filtering.""")

    def __init__(self, **params):
        super().__init__(**params)

        self.df = self.get_data()
        self.index = self.df.columns if self.axis else self.df.index
        self.names = self.names or self.index.names
        n_levels = self.n_levels or len(self.names)

        # either has a widgets dict or panel property
        self.widgets = {name: pn.widgets.Select(name=name) for name in self.names[:n_levels]}
        self.selectors = list(self.widgets.values())
        for selector in self.selectors:
            selector.param.watch(self._selector_changed, ['value'], onlychanged=True)

        options = list(self.index.get_level_values(0).unique())
        if self.empty_select:
            options = ['None'] + options
        self.selectors[0].options = options
        self.selectors[0].value = options[0]

    def update(self, *events):
        self.df = self.get_data()
        self.index = self.df.columns if self.axis else self.df.index  #assuming names stay the same
        super().update(*events)

    def _selector_changed(self, *events):
        for event in events:
            current_index = self.selectors.index(event.obj)  # Index of the selector which was changed

            try:  # try/except when we are at the last selector
                next_selector = self.selectors[current_index + 1]

                # Determine key/level to obtain new index which gives options for the next slider
                values = [selector.value for selector in self.selectors[:current_index + 1]]
                key = [value if value != 'None' else slice(None) for value in values]
                level = list(range(current_index+1))
                bools, current_columns = self.index.get_loc_level(key=key, level=level)
                options = list(current_columns.get_level_values(0).unique())
                if self.empty_select:
                    options = ['None'] + options
                next_selector.options = options
                if next_selector.value is None:  # If the selector was not set yet, set it to 'None'
                    next_selector.value = options[0]
            except IndexError:
                pass

        # set the df
        all_values = [selector.value for selector in self.selectors]
        self.key = tuple([value if value != 'None' else slice(None) for value in all_values])
        self.level = list(range(len(all_values)))

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        return dict(key=self.key, axis=self.axis, level=self.level, drop_level=self.drop_level)


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

