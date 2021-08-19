import logging
from itertools import groupby, count

import holoviews as hv
import pandas as pd
import panel as pn
import param
from bokeh.models import HoverTool, Span, Rect, Whisker
from bokeh.models.formatters import NumeralTickFormatter
from bokeh.plotting import figure
from lumen.filters import ParamFilter
from lumen.views import hvPlotView, View

from pyhdx.web.base import BokehFigurePanel, FigurePanel, MIN_BORDER_LEFT
from pyhdx.web.widgets import LoggingMarkdown, NGL


class hvPlotAppView(hvPlotView):

    def get_data(self):
        # get data filter using pandas query syntax?
        try:
            data = super().get_data()
        except (KeyError, ValueError) as e:
            #print(f'Empty data in {self.__class__}: {e}')
            return self.empty_df

        #data = super().get_data()

        if data.size > 2:  #todo fix this weird hack
            return data
        else:
            #print(f'got data but too small in  {self.__class__}, ')
            return self.empty_df

    @property
    def empty_df(self):
        dic = {self.x: [], self.y: []}
        if 'c' in self.kwargs:
            dic[self.kwargs['c']] = []
        return pd.DataFrame(dic)


class hvRectangleAppView(View):

    opts = param.Dict(default={}, doc="HoloViews option to apply on the plot.")

    view_type = 'rectangles'

    streaming = param.Boolean(default=False, doc="""
        Whether to stream new data to the plot or rerender the plot.""")

    def __init__(self, **params):
        # import hvplot.pandas # noqa
        # if 'dask' in sys.modules:
        #     try:
        #         import hvplot.dask # noqa
        #     except Exception:
        #         pass
        self._stream = None
        self._linked_objs = []
        super().__init__(**params)

    def get_panel(self):
        kwargs = self._get_params()
        #interactive? https://github.com/holoviz/panel/issues/1824
        return pn.pane.HoloViews(**kwargs)

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
        # processed = {}
        # for k, v in self.kwargs.items():
        #     if k.endswith('formatter') and isinstance(v, str) and '%' not in v:
        #         v = NumeralTickFormatter(format=v)
        #     processed[k] = v
        # if self.streaming:
        #     processed['stream'] = self._stream

        #hvplots stream? https://holoviews.org/user_guide/Streaming_Data.html

#        plot = hv.Rectangles([(0, 0, 1, 1), (2, 3, 4, 6), (0.5, 2, 1.5, 4), (2, 1, 3.5, 2.5)])

        processed = {}
        for k, v in self.kwargs.items():
            if k.endswith('formatter') and isinstance(v, str) and '%' not in v:
                v = NumeralTickFormatter(format=v)
            processed[k] = v
        if self.streaming:
            #processed['stream'] = self._stream

            plot = hv.DynamicMap(hv.Rectangles, streams=[self._stream])
            plot = plot.apply.opts(**self.opts) if self.opts else plot
        else:
            plot = hv.Rectangles(df)
            plot.opts(**self.opts) if self.opts else plot

        if self.selection_group or 'selection_expr' in self._param_watchers:
            plot = self._link_plot(plot)

        return plot

    def _get_params(self):
        df = self.get_data()
        if df is None:
            df = self.empty_df

        if self.streaming:
            from holoviews.streams import Pipe
            self._stream = Pipe(data=df)
        return dict(object=self.get_plot(df), sizing_mode='stretch_both')  # todo update sizing mode

    def get_data(self):
        #todo uniformify this method for all views
        try:
            return super().get_data()
        except (KeyError, ValueError) as e:
            #print(f'Empty data in {self.__class__}: {e}')
            return None
            #return self.empty_df

    def update(self, *events, invalidate_cache=True):
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
        own_parameters = [self.param[p] for p in self.param]
        own_events = events and all(
            isinstance(e.obj, ParamFilter) and
            (e.obj.parameter in own_parameters or
            e.new is self._ls.selection_expr)
            for e in events
        )
        if own_events:
            return False
        if invalidate_cache:
            self._cache = None
        if not self.streaming or self._stream is None:
            upd = self._update_panel()
            return upd
        if self.get_data() is not None:
            self._stream.send(self.get_data())
        return False


    @property
    def empty_df(self):
        return pd.DataFrame([[0] * 5], columns=['x0', 'x1', 'y0', 'y1', 'value'])


class NGLView(View):
    view_type = 'protein'

    #spin = param.Boolean(default=False)

    #js_files = {'ngl': "https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.37/dist/ngl.js"}

    def __init__(self, **params):
        super(NGLView, self).__init__(**params)
        self.ngl_view = NGL(sizing_mode='stretch_both')

    def get_panel(self):
        return self.ngl_view

    def update(self, *events, invalidate_cache=True):
        if invalidate_cache:
            self._cache = None
        data = self.get_data()
        if len(data.columns) > 1 or data.size < 1:
            # invalid number of columns
            self.ngl_view.color_list = [['white', "*"]]
        else:
            pd_series = data.iloc[:, 0]
            grp = pd_series.groupby(pd_series)

            color_list = []
            for c, pd_series in grp:
                result = [list(g) for _, g in groupby(pd_series.index, key=lambda n, c=count(): n - next(c))]

                resi = ' or '.join([f'{g[0]}-{g[-1]}' for g in result])
                color_list.append([c, resi])

            self.ngl_view.color_list = color_list

        # update panel?
        return self._update_panel()


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


class CoverageFigure(BokehFigurePanel):
    title = 'Coverage'
    accepted_tags = ['coverage']

    def __init__(self, *args, **params):
        super(CoverageFigure, self).__init__(*args, **params)

    def draw_figure(self):
        fig = figure(title=None, min_border=0, tools='pan,wheel_zoom,box_zoom,save,reset')
        fig.min_border_left = MIN_BORDER_LEFT
        fig.xaxis.axis_label = 'Residue number'

        return fig

    def render_sources(self, src_dict):  #todo , **render_kwargs
        tooltips = [('Pos', '$x{int}'),
                    ('Index', '@index'),
                    ('Start', '@start (@_start)'),
                    ('End', '@end (@_end)'),
                    ('Sequence', '@sequence'),
                    ('Score', '@scores'),
                    ('Uptake', '@uptake (@uptake_corrected / @ex_residues, @maxuptake)')]

        for name, data_source in src_dict.items():
            glyph = Rect(x='x', y='y', width='width', height=1, fill_color='color')
            renderer = self.figure.add_glyph(data_source.source, glyph)
            self.renderers[name] = renderer

            hovertool = HoverTool(renderers=[renderer], tooltips=tooltips)
            self.figure.add_tools(hovertool)


class LinearLogFigure(BokehFigurePanel):
    """base class for bokeh figure which can switch between log and linear axis
        This is a very slim base class (perhaps it should be a mixin?) (yes probably)
        and it should have a different name
    """

    def _redraw_event(self, event):
        """Redraw the figure with new kwargs passed to figure"""
        kwargs = {event.name: event.new.lower()}
        self.redraw(**kwargs)


class ThdFigure(LinearLogFigure):
    """base class for figures extending LinearLogFigure with threshold lines
    parent must have ClassificationControl
    """
    accepted_tags = ['mapping']

    # def __init__(self, parent, *args, **params):
    #     super(ThdFigure, self).__init__(parent, *args, **params)

    def draw_figure(self, **kwargs):
        y_axis_type = kwargs.pop('y_axis_type', 'log')
        x_axis_type = kwargs.pop('x_axis_type', 'linear')
        fig = figure(y_axis_type=y_axis_type, x_axis_type=x_axis_type,
                     tools='pan,wheel_zoom,box_zoom,save,reset', **kwargs)
        fig.min_border_left = MIN_BORDER_LEFT
        fig.xaxis.axis_label = 'Residue number'
        fig.yaxis.axis_label = self.y_label

        # todo refactor controller access to imporove upon current try/except and allow more controller agnostic behaviour
        try:
            for _ in range(self.control_panels['ClassificationControl'].param['num_colors'].bounds[1] - 1):
                sp = Span(location=0, dimension='width')
                sp.tags = ['thd']
                sp.visible = False
                fig.add_layout(sp)
        except KeyError:
            pass

        return fig

    def render_sources(self, src_dict, **render_kwargs):
        for name, data_source in src_dict.items():
            kwargs = {**data_source.render_kwargs, **render_kwargs}

            glyph_func = getattr(self.figure, data_source.renderer)
            renderer = glyph_func(source=data_source.source, name=name,
                                  legend_label=name, **kwargs)  #todo size is being specified at two different places now

            self.renderers[name] = renderer
            hovertool = HoverTool(renderers=[renderer],
                                  tooltips=[('Residue', '@r_number{int}'), (self.y_label, f'@{kwargs["y"]}')],
                                  mode='vline')
            self.figure.add_tools(hovertool)

        if self.renderers:
            self.figure.legend.click_policy = 'hide'

    def _draw_thds(self, *events):
        #todo duplicate code, subclass
        spans = self.figure.select(tags='thd')

        if not self.figure.renderers:
            return

        y_names, source_names = zip(*[(renderer.glyph.y, renderer.data_source.name) for renderer in self.figure.renderers])
        if not self.control_panels['ClassificationControl'].target in source_names:
            self._hide_thds()
            return
        if not self.control_panels['ClassificationControl'].quantity in y_names:
            self._hide_thds()
            return

        spans.sort(key=lambda x: x.id)
        for i, span in enumerate(spans):
            if i < len(self.control_panels['ClassificationControl'].values):
                span.location = self.control_panels['ClassificationControl'].values[i]
                span.visible = self.control_panels['ClassificationControl'].show_thds
            else:
                span.visible = False

    def _hide_thds(self):
        spans = self.figure.select(tags='thd')
        for i, span in enumerate(spans):
            span.visible = False


class RateFigure(ThdFigure):
    title = 'Rates'
    accepted_tags = [('mapping', 'rate')]
    y_label = 'Rate (min⁻¹)'

    def setup_hooks(self):
        super().setup_hooks()
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict, **render_kwargs)


class PFactFigure(ThdFigure):
    title = 'Protection Factors'
    accepted_tags = [('mapping', 'pfact')]
    y_label = 'Protection factor'

    def setup_hooks(self):
        #todo move to superclass?
        super().setup_hooks()
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict, y='pfact', **render_kwargs)


class DeltaGFigure(ThdFigure):
    title = 'DeltaG'
    accepted_tags = [('mapping', 'deltaG')]
    y_label = 'DeltaG (J/mol)'

    def setup_hooks(self):
        #todo move to superclass?
        super().setup_hooks()
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])

    def draw_figure(self, **kwargs):
        return super().draw_figure(y_axis_type='linear')

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict, y='deltaG', **render_kwargs)

        #todo make sure that if a new deltaG get plotted the graph is redrawn
        for name, data_source in src_dict.items():
            if 'covariance' in data_source.source.data.keys():
                # y = data_source.source.data['deltaG']
                # cov = data_source.source.data['covariance']
                # data_source.source.data['__upper'] = y + cov
                # data_source.source.data['__lower'] = y - cov

                whiskers = Whisker(source=data_source.source, base='r_number', upper='__upper', lower='__lower')
                self.figure.add_layout(whiskers)

    # def _data_updated_callback(self, attr, old, new):
    #     #get layouts and update the cds
    # #    print(attr, old, new)  # data, dict, propertyvalue columndata
    #
    #
    #     if not np.allclose(old['deltaG'], new['deltaG']):
    #         y = new['deltaG']
    #         # x = data_source.source.data['r_number']
    #         cov = new['covariance']
    #
    #         new['__upper'] = y + cov
    #         new['__lower'] = y - cov
    #
    #         super()._data_updated_callback(attr, old, new)


class BinaryComparisonFigure(ThdFigure):
    title = 'Binary Comparison'
    accepted_tags = [('comparison', 'mapping')]  # [ ('x' AND 'y') OR 'z' OR 'asdf']
    x_label = 'Residue number'  # move these to _redraw_kwargs?
    y_label = 'Difference'

    def __init__(self, parent, *args, **params):
        super(BinaryComparisonFigure, self).__init__(parent, *args, **params)
        self._redraw_kwargs['y_axis_type'] = 'linear'

    @property
    def y_kwarg(self):
        try:
            plot_quantity = self.control_panels['ClassificationControl'].quantity
            if plot_quantity is not None:
                return plot_quantity
        except KeyError:
            pass
        return None

    def render_sources(self, src_dict, **render_kwargs):
        kwargs = {**self._render_kwargs, **render_kwargs}
        super().render_sources(src_dict, **kwargs)

    def setup_hooks(self):
        super().setup_hooks()
        #todo this should be resolved in some other way than the name
        try:
            self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])
            self.control_panels['ClassificationControl'].param.watch(self._log_space_updated, ['log_space'])
            self.control_panels['ClassificationControl'].param.watch(self._quantity_updated, ['quantity'])
        except KeyError:
            pass

        try:
            self.control_panels['ColoringControl'].param.watch(self._draw_thds, ['values', 'show_thds'])
            self.control_panels['ColoringControl'].param.watch(self._log_space_updated, ['log_space'])
        except KeyError:
            pass

    #todo group into one function?
    def _quantity_updated(self, event):
        self._render_kwargs['y'] = event.new
        self.y_label = event.new
        self.redraw(**self._redraw_kwargs) # todo auto pass kwargs?

    # def _target_updated(self, event):  #event. cls, name, new, obj, old, type, what
    #     #todo make redraw accept events
    #     self._render_kwargs['target'] = event.new
    #     self.redraw(**self._redraw_kwargs) # todo auto pass kwargs?

    def _log_space_updated(self, event):
        if event.new:   # True, log space
            self._redraw_kwargs['y_axis_type'] = 'log'
        else:
            self._redraw_kwargs['y_axis_type'] = 'linear'
        self.redraw(**self._redraw_kwargs) # todo auto pass kwargs?

    def draw_figure(self, **kwargs):
        y_axis_type = kwargs.pop('y_axis_type', 'linear')
        return super().draw_figure(y_axis_type=y_axis_type, **kwargs)


class ScoresFigure(LinearLogFigure):
    title = 'Scores'
    accepted_tags = [('scores', 'mapping')]
    x_label = 'Residue number'
    y_label = 'Value'

    def render_sources(self, src_dict):
        for name, data_source in src_dict.items():
            for y_field in data_source.render_kwargs['y']:
                glyph_func = getattr(self.figure, data_source.renderer)
                kwargs = data_source.render_kwargs.copy()
                kwargs.pop('y')

                renderer = glyph_func(**kwargs, y=y_field, source=data_source.source, name=name,
                                      legend_label=f'{y_field}')

                self.renderers[name] = renderer
                # hovertool = HoverTool(renderers=[renderer],
                #                       tooltips=[('Residue', '@r_number{int}'), (self.y_label, f'@{data_source.render_kwargs["y"]}')],
                #                       mode='vline')
                # self.figure.add_tools(hovertool)

            if self.renderers:
                self.figure.legend.click_policy = 'hide'


class SingleValueFigure(LinearLogFigure):
    title = 'Values'
    accepted_tags = [('comparison', 'mapping')]
    x_label = 'Residue number'
    y_label = 'Value'

    def render_sources(self, src_dict):
        for name, data_source in src_dict.items():
            for field, render_func in zip(['value1', 'value2'], ['triangle', 'square']):
                glyph_func = getattr(self.figure, render_func)
                kwargs = data_source.render_kwargs.copy()
                kwargs.pop('y')

                renderer = glyph_func(**kwargs, y=field, source=data_source.source, name=name,
                                      legend_label=name + f'_{field}')

                self.renderers[name] = renderer
                hovertool = HoverTool(renderers=[renderer],
                                      tooltips=[('Residue', '@r_number{int}'), (self.y_label, f'@{data_source.render_kwargs["y"]}')],
                                      mode='vline')
                self.figure.add_tools(hovertool)

            if self.renderers:
                self.figure.legend.click_policy = 'hide'

    def draw_figure(self, **kwargs):
        y_axis_type = kwargs.pop('y_axis_type', 'linear')
        return super().draw_figure(y_axis_type=y_axis_type, **kwargs)


class SingleFigure(BinaryComparisonFigure):
    title = 'Values'
    y_label = 'Value'


class FitResultFigure(BokehFigurePanel):
    title = 'Fit Result'
    accepted_tags = ['uptake_curve']
    y_label = 'Uptake corrected'
    x_label = 'Time'

    def __init__(self, parent, *args, **params):
        super(FitResultFigure, self).__init__(parent, *args, **params)

    def setup_hooks(self):
        self.control_panels['FitResultControl'].param.watch(self._redraw_event, ['x_axis_type'])

    def _redraw_event(self, event):
        """Redraw the figure with new kwargs passed to figure"""
        kwargs = {event.name: event.new.lower()}
        self.redraw(**kwargs)

    def draw_figure(self, **kwargs):
        fig = super().draw_figure(x_axis_type=self.control_panels['FitResultControl'].x_axis_type.lower())
        return fig

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict)

        if self.renderers:
            self.figure.legend.location = "bottom_right"


class LoggingFigure(FigurePanel):
    title = 'Logging'

    def __init__(self, *args, **params):
        super(LoggingFigure, self).__init__(*args, **params)
        self.markdown = LoggingMarkdown('### Log Window \n', sizing_mode='stretch_both')

        sh = logging.StreamHandler(self.markdown)
        sh.terminator = '  \n'
        formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', "%Y-%m-%d %H:%M:%S")
        sh.setFormatter(formatter)
        #sh.setLevel(logging.DEBUG)
        self.parent.logger.addHandler(sh)

    def setup_hooks(self):
        # todo add function on base class for doing these things (with try/except) ?
        # Or the users has to be responsible and set this up correctly
        # if the hook is unwanted, the class should be subclassed with override on setup_hooks

        try:
            self.parent.control_panels['OptionsControl'].param.watch(self._update_log_level, ['log_level'])
            self.parent.control_panels['OptionsControl'].param.trigger('log_level')
        except KeyError:
            self.parent.logger.debug('Control panel OptionsControl not founc')

    def _update_log_level(self, event):
        self.parent.logger.setLevel(event.new)

    @property
    def panel(self):
        return self.markdown


class ImageFigure(LinearLogFigure):
    title = 'Image'

    accepted_tags = ['image']
    x_label = 'Residue number'
    y_label = 'Time (probably)'

    def draw_figure(self, **kwargs):
        figure = super().draw_figure()
        figure.x_range.range_padding = figure.y_range.range_padding = 0

        return figure

    def render_sources(self, src_dict, **render_kwargs):
        render_kwargs.pop('color', None)
        for name, data_source in src_dict.items():
            ds_kwargs = data_source.render_kwargs
            ds_kwargs.pop('color', None)   # todo fix color winding up here int he first place
            renderer = self.figure.image_rgba(source=data_source.source, image='img', **ds_kwargs)
            hovertool = HoverTool(renderers=[renderer],
                                  tooltips=[("x", "$x"), ("y", "$y"), ("value", "@scores")]
                                  )
            self.figure.add_tools(hovertool)
    #
