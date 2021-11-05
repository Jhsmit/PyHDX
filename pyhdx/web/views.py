import itertools
import logging
from functools import partial
from itertools import groupby, count

import holoviews as hv
from holoviews.streams import Pipe
import pandas as pd
import panel as pn
import param
from bokeh.models import HoverTool, Span, Rect, Whisker
from bokeh.models.formatters import NumeralTickFormatter
from bokeh.plotting import figure
from lumen.filters import ParamFilter
from lumen.views import hvPlotView, View
from panel.pane.base import PaneBase

from pyhdx.support import autowrap
from pyhdx.web.base import BokehFigurePanel, FigurePanel, MIN_BORDER_LEFT
from pyhdx.web.filters import AppFilter
from pyhdx.web.sources import AppSource
from pyhdx.web.widgets import LoggingMarkdown, NGL, REPRESENTATIONS, COLOR_SCHEMES

import numpy as np


class AppView(param.Parameterized):
    """Base view object.

    Inspired by Holoviz Lumen's View objects"""

    # filters = param.List(constant=True, doc="""
    #     A list of Filter object providing the query parameters for the
    #     Source.""")

    source = param.ClassSelector(class_=(AppSource, AppFilter),
                                 constant=True,
                                 precedence=-1,
                                 doc="""
        The Source to query for the data.""")

    opts = param.Dict(default={}, doc="HoloViews option to apply on the plot.",
                      precedence=-1,
                      constant=True)

    dependencies = param.List(
        default=[],
        precedence=-1,
        doc="Additional dependencies which trigger update when their `updated` event fires")

    view_type = None

    def __init__(self, **params):
        super().__init__(**params)
        # todo allow for kwargs to be passed to DynamicMap's func

        for dep in self.dependencies:
            dep.param.watch(self._update, ['updated'])

        self.widgets = self.make_dict()

        self._panel = None
        self._updates = None  # what does this do?

    def _update(self, *events):
        self.update()  # todo or just catch events in the update function?  (probably this)

    def make_dict(self):
        #todo baseclass filter/controller/view with the same widget generation logic?
        """dict of widgets to be shown
        override this method to get custom mapping

        """
        return self.generate_widgets()

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""

        names = [p for p in self.param if self.param[p].precedence is None or self.param[p].precedence > 1]
        widgets = pn.Param(self.param, show_name=False, show_labels=True, widgets=kwargs)

        return {k: v for k, v in zip(names[1:], widgets)}

    def get_data(self):
        """
        Queries the Source for the specified table applying any
        filters and transformations specified on the View. Unlike
        `get_value` this should be used when multiple return values
        are expected.

        Returns
        -------
        DataFrame
            The queried table after filtering and transformations are
            applied.
        """
        # if self._cache is not None:
        #     return self._cache

        df = self.source.get()

        return df

    def _update_panel(self, *events):
        """
        Updates the cached Panel object and returns a boolean value
        indicating whether a rerender is required.
        """
        if self._panel is not None:
            self._cleanup()
            self._updates = self._get_params()
            if self._updates is not None:
                return False
        self._panel = self.get_panel()
        return True


class hvAppView(AppView):

    def __init__(self, **params):
        super().__init__(**params)
        self._stream = None

    @param.depends('source.updated', watch=True)  # no watch? # todo cache / checking if updates are needed?
    def update(self):
        """
        Triggers an update in the View.

        Parameters
        ----------
        events: tuple
            param events that may trigger an update.
        invalidate_cache : bool
            Whether to clear the View's cache.

        Returns
        -------
        stale : bool
            Whether the panel on the View is stale and needs to be
            rerendered.
        """
        # Skip events triggered by a parameter change on this View
        # own_parameters = [self.param[p] for p in self.param]
        # own_events = events and all(
        #     isinstance(e.obj, ParamFilter) and
        #     (e.obj.parameter in own_parameters or
        #     e.new is self._ls.selection_expr)
        #     for e in events
        # )
        # if own_events:
        #     return False
        # if invalidate_cache:
        #     self._cache = None
        if self._stream is None:
            return self._update_panel()
        if self.get_data() is not None:
            self._stream.send(self.get_data())
        return False

    def get_panel(self):
        kwargs = self._get_params()
        #interactive? https://github.com/holoviz/panel/issues/1824
        return pn.pane.HoloViews(**kwargs)

    def _get_params(self):
        df = self.get_data()

        self._stream = Pipe(data=df)
        return dict(object=self.get_plot(df), sizing_mode='stretch_both')  # todo update sizing mode

    @property
    def panel(self):
        if isinstance(self._panel, PaneBase):
            pane = self._panel
            if len(pane.layout) == 1 and pane._unpack:
                return pane.layout[0]
            return pane._layout
        return self._panel


class hvScatterAppView(hvAppView):

    view_type = 'scatter'

    x = param.String(doc="The column to render on the x-axis.")  # todo these should be selectors

    y = param.String(doc="The column to render on the y-axis.")

    def __init__(self, **params):
        self._stream = None
        super().__init__(**params)

    def get_plot(self, df):
        """

        Parameters
        ----------
        df

        Returns
        -------

        """

        func = partial(hv.Scatter, kdims=[self.x], vdims=[self.y])
        plot = hv.DynamicMap(func, streams=[self._stream])
        plot = plot.apply.opts(**self.opts) if self.opts else plot

        return plot


    # @property
    # def empty_df(self):
    #     dic = {self.x: [], self.y: []}
    #     if 'c' in self.kwargs:
    #         dic[self.kwargs['c']] = []
    #     return pd.DataFrame(dic)


class hvRectangleAppView(hvAppView):

    view_type = 'rectangles'

    left = param.String('start', doc="Field name to use for left coordinate")

    right = param.String('stop', doc="Field name to use for for the right coordinate")

    height = param.Integer(1, doc="Height of the rectangles", constant=True)

    passthrough = param.List()

    wrap = param.Integer(None, doc="Number of peptides to go down on y axis before wrapping back around to the top")
    step = param.Integer(5, bounds=(1, None), doc="Step size used for finding 'wrap' when its not specified")
    margin = param.Integer(4, doc="Margin space to keep between peptides when finding 'wrap'")

    def __init__(self, **params):

        #todo left and right cannot be none?
        super().__init__(**params)

    def get_data(self):
        df = super().get_data()

        #todo this nonsense below should be the action some kind of filter

        wrap = self.wrap or autowrap(df[self.left].to_numpy(dtype=int), df[self.right].to_numpy(dtype=int), margin=self.margin, step=self.step)

        # Order of columns determines their role, not the names
        columns = ['x0', 'y0', 'x1', 'y1']  # bottom-left (x0, y0) and top right (x1, y1)
        output_table = pd.DataFrame(index=df.index, columns=columns)
        output_table['x0'] = df[self.left] - 0.5
        output_table['x1'] = df[self.right] - 0.5
        cycle = itertools.cycle(range(self.height*wrap, 0, -self.height))  # Starts at y top, cycles with wrap
        yvals = np.array(list(itertools.islice(cycle, len(df)))) # Repeat cycle until the correct length
        output_table['y0'] = yvals - self.height
        output_table['y1'] = yvals

        if self.passthrough:
            for item in self.passthrough:
                assert item not in ['value', 'index'], "Invalid field name, 'index' and 'value' names are reserved"
                output_table[item] = df[item]
        output_table['index'] = df.index

        return output_table

    def get_plot(self, df):
        """
        Dataframe df must have columns x0, y0, x1, y1 (in this order) for coordinates
        bottom-left (x0, y0) and top right (x1, y1). Optionally a fifth value-column can be provided for colors

        Parameters
        ----------
        df

        Returns
        -------

        """

        plot = hv.DynamicMap(hv.Rectangles, streams=[self._stream])
        plot = plot.apply.opts(**self.opts) if self.opts else plot

        return plot

    @property
    def empty_df(self):
        return pd.DataFrame([[0] * 5], columns=['x0', 'x1', 'y0', 'y1', 'value'])


class NGLView(AppView):
    view_type = 'ngl'

    # todo additioal render options (fog etc)

    representation = param.Selector(default='cartoon', objects=REPRESENTATIONS)

    spin = param.Boolean(False)

    color_scheme = param.Selector(default='custom', objects=COLOR_SCHEMES)

    custom_color_scheme = param.List(precedence=-1)

    background_color = param.Color(default='#F7F7F7')

    object = param.String('', doc='pdb string object', precedence=-1)

    def __init__(self, **params):
        super(NGLView, self).__init__(**params)
        self.ngl_view = NGL(sizing_mode='stretch_both', #todo sanitize order
                            extension='pdb',
                            spin=self.spin,
                            color_scheme=self.color_scheme,
                            representation=self.representation,
                            object=self.object,

                            )

        params = self.param.params().keys() & self.ngl_view.param.params().keys() - {'name'}
        self.param.watch(self._update_params, list(params))
        # todo is there a better way to couple two params?

    def _update_params(self, *events):
        for event in events:
            setattr(self.ngl_view, event.name, event.new)

    def get_panel(self):
        return self.ngl_view

    def _cleanup(self):
        return None

    def _get_params(self):
        return None

    @param.depends('source.updated', watch=True)
    def update(self):
        colors = self.get_data()
        if colors is not None:
            grp = colors.groupby(colors)

            color_list = []
            for c, pd_series in grp:
                result = [list(g) for _, g in groupby(pd_series.index, key=lambda n, c=count(): n - next(c))]

                resi = ' or '.join([f'{g[0]}-{g[-1]}' for g in result])
                color_list.append([c, resi])

            self.custom_color_scheme = color_list

        return self._update_panel()  # what does this do? redraws if not redrawn?

    @property
    def panel(self):  # why the panebase unpack?
        if isinstance(self._panel, PaneBase):
            pane = self._panel
            if len(pane.layout) == 1 and pane._unpack:
                return pane.layout[0]
            return pane._layout
        return self._panel


class LoggingView(View):
    view_type = 'logging'

    logger = param.ClassSelector(logging.Logger, doc='Logger object to show in Log view')

    level = param.Integer(default=10, doc='Logging level of the streamhandler redirecting logs to the view')

    def __init__(self, *args, **params):
        super(LoggingView, self).__init__(**params)
        self.markdown = LoggingMarkdown('### Log Window \n', sizing_mode='stretch_both')

        self.sh = logging.StreamHandler(self.markdown)
        self.sh.terminator = '  \n'
        self.sh.setLevel(self.level)
        formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', "%Y-%m-%d %H:%M:%S")
        self.sh.setFormatter(formatter)
        self.logger.addHandler(self.sh)

    @param.depends('level', watch=True)
    def _level_updated(self):
        self.sh.setLevel(self.level)

    @property
    def panel(self):
        return self.markdown

    def update(self, *events, invalidate_cache=True):
        pass
