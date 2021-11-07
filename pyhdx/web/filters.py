import itertools

import numpy as np
import pandas as pd
import panel as pn
import param
from param.parameterized import default_label_formatter

from pyhdx.support import autowrap
from pyhdx.web.sources import AppSourceBase


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


class TableSourceFilter(AppFilterBase):  #todo rename to something that includes table
    """filter which picks the correct table from the source"""

    _type = 'table_source'

    source = param.ClassSelector(class_=AppSourceBase)

    table = param.Selector(default=None, doc="""
      The table being filtered. """)

    def __init__(self, table_options=None, **params):
        self.table_options = table_options
        super().__init__(**params)

        self.widgets = {'table': pn.pane.panel(self.param.table)}

    #todo allow auto generate widgets as in control panels /  views

    def get(self):
        df = self.source.get(self.table)  # returns None on KeyError #todo change to source.get_table
        return df

    @param.depends('table', watch=True)
    def _table_updated(self):
        self.updated = True

    @param.depends('source.updated', watch=True)  # todo split
    def update(self):
        options = self.source.get_tables()
        if self.table_options:
            options = [t for t in options if t in self.table_options]
        #with param.parameterized.discard_events(self):  # no triggers from setting table, manual with updated=True
        self.param['table'].objects = options
        if not self.table and options:  # todo more nonsense as in filters?
            self.table = options[0]
        self.updated = True


class AppFilter(AppFilterBase):
    """filter which acts on previous filters in a chain. source is also afilter"""
    source = param.ClassSelector(class_=AppFilterBase)

    def get(self):
        df = self.source.get()

        return df

    @param.depends('source.updated', watch=True)
    def update(self):
        self.updated = True


class CrossSectionFilter(AppFilter):

    _type = 'cross_section'

    pd_function = 'xs'

    key = param.List()

    axis = param.Integer(1, bounds=[0, 1])

    n_levels = param.Integer(0, doc="Number of levels. negative to count from the back of the list")

    level = param.List()

    drop_level = param.Boolean(True)

    names = param.List(None, doc="List of label names for widgets")

    empty_select = param.Boolean(default=False, doc="""
        Add an option to Select widgets to indicate no filtering.""")

    def __init__(self, **params):
        super().__init__(**params)
        self.index = None  # index is the df index which determines the selector's options
        self.update()

    @param.depends('source.updated', watch=True)
    def update(self):
        #todo only redraw if only options are changed or always?
        #todo remove watchers when new filters are created?

        old_index = self.index
        df = self.source.get()
        if df is None:
            return

        self.index = df.columns if self.axis else df.index
        self._names = self.names or self.index.names

        if old_index is not None and self.index.nlevels == old_index.nlevels:
            # no redraw needed, only update selectors options
            options = list(self.index.unique(level=0))
            self.selectors[0].options = options
            self.selectors[0].param.trigger('value')
            for name, selector in zip(self._names, self.selectors):
                selector.name = name  # todo requires testing if the names are really updated or not (they arent)
                selector.label = name  # todo requires testing if the names are really updated or not
                self.redrawn = True
        else:
            self.redraw()

        self.updated = True

    def redraw(self):
        # create new widgets
        if self.n_levels <= 0:
            n_levels = len(self._names) + self.n_levels
        else:
            n_levels = self.n_levels

        self.widgets = {name: pn.widgets.Select(name=default_label_formatter(name)) for name in self._names[:n_levels]}

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
        if df is None:
            return df
        else:
            kwargs = self.pd_kwargs
            df = df.xs(**self.pd_kwargs)
            return df

    def _selector_changed(self, *events):
        #this sends multiple updated events as it triggers changes in other selectors
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
        self.key = [value if value != 'None' else slice(None) for value in all_values]
        self.level = list(range(len(all_values)))

        self.updated = True

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        return dict(key=tuple(self.key), axis=self.axis, level=self.level, drop_level=self.drop_level)


class ApplyCmapOptFilter(AppFilter):

    _type = 'apply_cmap_opt'

    opts = param.Selector(doc='cmap opts list to choose from',
                          label='Color transform', objects=[],
                          )  #todo refactor to cmap_opt? color transform?

    #@clasmethod
    #def check_args(... )  #todo method for constructor to see if the supplied kwargs are correct for this object

    def __init__(self, opts, **params): #opts: list of opts objects
        self._opts_dict = {o.name: o for o in opts}
        opts = list(self._opts_dict.keys())
        params['opts'] = opts
        super().__init__(**params)
        self.widgets = {'opts': pn.pane.panel(self.param.opts)}

    def get(self):
        df = self.source.get()  #todo refactor df to data as it can also be a series?
        if df is None:
            return None

        if df.columns.size != 1:
            raise ValueError("Invalid number of columns, must be 1")

        pd_series = df[df.columns[0]]
        opts_obj = self._opts_dict[self.opts]
        if pd.api.types.is_numeric_dtype(pd_series):
            colors = opts_obj.apply(pd_series)

            return colors
        else:
            return None

    @param.depends('opts', watch=True)
    def _opts_changed(self):
        self.updated = True  # opts options are the same but selection changed, signal

    @param.depends('source.updated', watch=True)
    def update(self):
        pd_series = self.source.get()
        #todo just show all, later deal with setting the correct one? (infer from previous filter setting)
        if pd_series is None:
            # with param.parameterized.discard_events(self):
            self.param['opts'].objects = []
            self.opts = None
        else:
            options = list(self._opts_dict.keys())  # this needs updating as opts_dict is static
            # with turn off param triggers, then update (unexpected results)
           # with param.parameterized.discard_events(self):
            self.param['opts'].objects = options  # TODO: pr/issue:? when setting objects which does not include the current setting selector is not reset?
            if self.opts is None and options:
                self.opts = options[0]
            elif self.opts not in options and options:  #todo or
                self.opts = options[0]
            elif not options:
                self.opts = None

        self.updated = True


class RectangleLayoutFilter(AppFilter):
    """
    Takes data with peptide start, end positions and transforms it to bottom, left, right, top positions of rectangles.

    By default, for a given interval start, end (inclusive, exclusive) a rectangle is returned with a height of spanning from
    (start - 0.5) to (end - 0.5)

    """

    _type = 'rectangle_layout'

    left = param.String('start', doc="Field name to use for left coordinate")
    right = param.String('end', doc="Field name to use for for the right coordinate")

    height = param.Integer(1, doc="Height of the rectangles", constant=True)
    value = param.String(None, doc="Optional field name to pass through as value column")
    passthrough = param.List([], doc="Optional field names to pass through (for hovertools)")

    wrap = param.Integer(doc="Amount of peptides to plot on the y axis before wrapping around to y axis top")
    step = param.Integer(5, bounds=(1, None), doc="Step size used for finding 'wrap' when its not specified")
    margin = param.Integer(4, doc="Margin space to keep between peptides when finding 'wrap'")

    def get(self):
        df = self.source.get()
        if df is None:
            return None

        if not self.wrap:
            wrap = autowrap(df[self.left].to_numpy(dtype=int), df[self.right].to_numpy(dtype=int),  # todo fix types to be ints in the df
                            margin=self.margin, step=self.step)

        # Order of columns determines their role, not the names
        columns = ['x0', 'y0', 'x1', 'y1']  # bottom-left (x0, y0) and top right (x1, y1)
        output_table = pd.DataFrame(index=df.index, columns=columns)
        output_table['x0'] = df[self.left] - 0.5
        output_table['x1'] = df[self.right] - 0.5
        cycle = itertools.cycle(range(self.height*wrap, 0, -self.height))  # Starts at y top, cyles with wrap
        yvals = np.array(list(itertools.islice(cycle, len(df)))) # Repeat cycle until the correct length
        output_table['y0'] = yvals - self.height
        output_table['y1'] = yvals


        if self.passthrough:
            for item in self.passthrough:
                assert item not in ['value', 'index'], "Invalid field name, 'index' and 'value' names are reserved"
                output_table[item] = df[item]
        output_table['index'] = df.index

        return output_table


class GenericFilter(AppFilter):
    _type = 'generic'

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


class RescaleFilter(GenericFilter):
    """Rescale a single column"""

    _type = 'rescale'

    pd_function = 'assign'

    column = param.String(doc='Name of the column to rescale')

    scale_factor = param.Number(1.)

    def get(self):  # todo perhaps some kind of decorator that returns nonealwasy?
        df = self.source.get()
        if df is None:
            return None
        df = df.assign(**self.pd_kwargs)

        return df

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        return {self.column: lambda x: x[self.column]*self.scale_factor}




class PivotFilter(GenericFilter):
    _type = 'pivot'

    pd_function = 'pivot'

    index = param.ClassSelector(class_=[str, list])

    columns = param.ClassSelector(class_=[str, list])

    values = param.ClassSelector(class_=[str, list])

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        #todo get_params func which finds the correct params here
        return dict(index=self.index, columns=self.columns, values=self.values)

class StackFilter(GenericFilter):
    _type = 'stack'

    pd_function = 'stack'

    level = param.Integer(-1)  #actually int, str, or list of these, default -1 (last level)
    #level = param.ClassSelector(_class=[int, str, list])  # where list is list of (str | int)

    dropna = param.Boolean(True)

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function to apply filter"""
        #todo get_params func which finds the correct params here
        return dict(level=self.level, dropna=self.dropna)



class ConcatFilter(AppFilter):
    """app filter which combines multiple dfs"""
    def __init__(self, **params):
        raise NotImplementedError()

class TransformFilter(AppFilter):
    pd_function = param.String('transform')
    def __init__(self, **params):
        raise NotImplementedError()