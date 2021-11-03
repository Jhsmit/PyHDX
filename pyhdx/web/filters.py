from lumen.filters import Filter
from lumen.sources import Source
from param.parameterized import default_label_formatter
import param

from pyhdx.web.sources import AppSource
from pyhdx.web.widgets import ColoredStaticText
import panel as pn

class AppFilterBase(param.Parameterized):

    """these filters get the data from source"""



    widgets = param.Dict(default={})

    updated = param.Event()

    redrawn = param.Event(doc="event gets triggered when widgets are changed and the controller needs to redraw them")

    def __init__(self, **params):
        super().__init__(**params)

    def get(self):
        """method called to get the dataframe"""
        return None


class AppSourceFilter(AppFilterBase):
    """filter which picks the correct table from the source"""

    source = param.ClassSelector(class_=AppSource)

    table = param.Selector(default=None, doc="""
      The table being filtered. """)

    def __init__(self, **params):
        super().__init__(**params)

        self.widgets = {'table': pn.pane.panel(self.param.table)}

    def get(self):
        df = self.source.get(self.table)
        return df

    @param.depends('source.updated', 'table', watch=True)
    def update(self):
        self.updated = True


class AppFilter(AppFilterBase):
    """filter which acts on previous filters in a chain. source is also afilter"""
    source = param.ClassSelector(class_=AppFilterBase)

    def get(self):
        df = self.source.get()

        return df


class AppFilterOld(AppFilterBase):
    """

    """

    widgets = param.Dict(default={})

    pd_function = param.String(doc='Pandas function which this filter applies to the DataFrame ')

    def __init__(self, table_select=False, **params):
        super().__init__(**params)

        if table_select:
            self._update_table_select()
            self.widgets = {'table': pn.pane.panel(self.param.table)}
            self.param.watch(self.update, 'table')

        # self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        # super().__init__(**{k: v for k, v in params.items() if k in self.param})

    def get_data(self):
        """Equivalent to view's get_data method"""

        queries = [filt.query for filt in self.filters]
        data = self.source.get(self.table, *queries)
        for transform in self.transforms:
            data = transform.apply(data)

        return data

    def _update_table_select(self):
        options = self.source.get_tables()
        self.param['table'].objects = options
        if not self.table:
            self.table = options[0]

    @param.depends('source.updated', watch=True)
    def update(self):
        #todo table spec
        """gets called when the source event triggers"""

        self._update_table_select()

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


class CrossSectionFilter(AppFilter):

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
        self.index = None  # index is the df index which determines the selector's options
        # or just call update on init?
        self.update()
        # df = self.get_data()
        # self.index = None#self.df.columns if self.axis else self.df.index
        # self._names = self.names or self.index.names
        #
        # # either has a widgets dict or panel property  ( probably we're going for widgets only?)
        # widgets = {name: pn.widgets.Select(name=name) for name in self.names[:n_levels]}
        # self.widgets.update(**widgets)
        # self.selectors = list(widgets.values())  # selectors is selectors minus the optional table selector
        # for selector in self.selectors:
        #     selector.param.watch(self._selector_changed, ['value'], onlychanged=True)
        #
        # options = list(self.index.get_level_values(0).unique())
        # if self.empty_select:
        #     options = ['None'] + options
        # self.selectors[0].options = options
        # self.selectors[0].value = options[0]

    #todo should rerender widgets upon source updated (perhaps)

    @param.depends('source.updated', watch=True)
    def update(self):
        #todo only redraw if only options are changed

        old_index = self.index
        df = self.source.get()
        self.index = df.columns if self.axis else df.index
        self._names = self.names or self.index.names

        if old_index is not None and self.index.nlevels == old_index:
            # no redraw needed, only update selectors options
            options = list(self.index.unique(level=0))
            self.selectors[0].options = options
            self.selectors[0].trigger('value')  # is this how it works?
            for name, selector in zip(self._names, self.selecotors):
                selector.name = name  # todo requires testing if the names are really updated or not
        else:
            self.redraw()

        self.updated = True

    def redraw(self):
        # create new widgets
        n_levels = self.n_levels or len(self._names)
        self.widgets = {name: pn.widgets.Select(name=name) for name in self._names[:n_levels]}
        #self.widgets.update(**widgets)
        self.selectors = list(self.widgets.values())
        for selector in self.selectors:
            selector.param.watch(self._selector_changed, ['value'], onlychanged=True)

        options = list(self.index.get_level_values(0).unique())
        if self.empty_select:  # todo use Nonetype? -> allow_none kwarg for Select?
            options = ['None'] + options
        self.selectors[0].options = options
        self.selectors[0].value = options[0]

        self.redrawn = True

    #todo cache df?
    def get(self):
        df = self.source.get()
        df = df.xs(**self.pd_kwargs)

        return df

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
                options = list(current_columns.get_level_values(0).unique())   # current_columns.unique(level=0)
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

        #signal the change
        self.updated = True

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        return dict(key=self.key, axis=self.axis, level=self.level, drop_level=self.drop_level)


class TransformFilter(AppFilter):
    pd_function = param.String('transform')


class StackFilter(AppFilter):
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

class GenericFilter(AppFilter):
    pd_function = param.String()

    def __init__(self, **params):
        self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**{k: v for k, v in params.items() if k in self.param})

    def get(self):
        df = self.source.get()
        func = getattr(df, self.pd_function)
        df = func(**self.pd_kwargs)

        return df

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas functin to apply filter"""
        return self.kwargs


class RescaleFilter(AppFilter):
    """Rescale a single column"""

    pd_function = 'assign'

    column = param.String(doc='Name of the column to rescale')

    scale_factor = param.Number(1.)

    def get(self):
        df = self.source.get()
        df = df.assign(**self.pd_kwargs)

        return df

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        return {self.column: lambda x: x[self.column]*self.scale_factor}


class PivotFilter(AppFilter):
    pd_function = 'pivot'

    index = param.ClassSelector(class_=[str, list])

    columns = param.ClassSelector(class_=[str, list])

    values = param.ClassSelector(class_=[str, list])

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        #todo get_params func which finds the correct params here
        return dict(index=self.index, columns=self.columns, values=self.values)


class AppWidgetFilter(AppFilter):

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


class UniqueValuesFilter(AppWidgetFilter):
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


class MultiIndexSelectFilter(AppWidgetFilter):
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

