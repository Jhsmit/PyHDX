from .base import BokehFigurePanel, FigurePanel, DEFAULT_RENDERERS, DEFAULT_COLORS, MIN_BORDER_LEFT
from .widgets import NGLViewer, LoggingMarkdown
from .log import setup_md_log
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column
from bokeh.models import LabelSet, ColumnDataSource, HoverTool, GlyphRenderer, Span, Rect, Range1d
from bokeh.models.markers import Triangle, Circle, Diamond
import panel as pn
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import logging

import param


class CoverageFigure(BokehFigurePanel):
    title = 'Coverage'
    accepted_sources = ['coverage']

    def __init__(self, *args, **params):
        super(CoverageFigure, self).__init__(*args, **params)

    def draw_figure(self):
        fig = figure(title=None, min_border=0, tools='pan,wheel_zoom,box_zoom,save,reset')
        fig.min_border_left = MIN_BORDER_LEFT
        fig.xaxis.axis_label = 'Residue number'

        return fig

    def render_sources(self, src_dict):
        tooltips = [('Pos', '$x{int}'),
                    ('Index', '@index'),
                    ('Start', '@start (@_start)'),
                    ('End', '@end (@_end)'),
                    ('Sequence', '@sequence'),
                    ('Score', '@scores'),
                    ('Uptake', '@uptake (@uptake_corrected / @ex_residues, @maxuptake)')]

        for name, source in src_dict.items():
            glyph = Rect(x='x', y='y', width='width', height=1, fill_color='color')
            renderer = self.figure.add_glyph(source, glyph)
            self.renderers[name] = renderer

            hovertool = HoverTool(renderers=[renderer], tooltips=tooltips)
            self.figure.add_tools(hovertool)


class ThdLogFigure(BokehFigurePanel):
    """base class for pfact / rates figure panels which both feature y log axis and thresholding"""

    def __init__(self, parent, *args, **params):
        super(ThdLogFigure, self).__init__(parent, *args, **params)
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])

    def draw_figure(self):
        fig = figure(y_axis_type='log', tools='pan,wheel_zoom,box_zoom,save,reset')
        fig.min_border_left = MIN_BORDER_LEFT
        fig.xaxis.axis_label = 'Residue number'
        fig.yaxis.axis_label = self.y_label

        for _ in range(self.control_panels['ClassificationControl'].param['num_colors'].bounds[1] - 1):  # todo refactor controller access
            sp = Span(location=0, dimension='width')
            sp.tags = ['thd']
            sp.visible = False
            fig.add_layout(sp)

        return fig

    def render_sources(self, src_dict):
        #todo perhaps refactor to single source
        for name, source in src_dict.items():
            func_name = DEFAULT_RENDERERS[name]
            glyph_func = getattr(self.figure, func_name)
            renderer = glyph_func(x='r_number', y='y', color='color', source=source, legend_label=name,
                                  size=10, name=name)
            self.renderers[name] = renderer
            hovertool = HoverTool(renderers=[renderer], tooltips=[('Residue', '@r_number{int}'), (self.y_label, '@y')],
                                  mode='vline')
            self.figure.add_tools(hovertool)

        if self.renderers:
            self.figure.legend.click_policy = 'hide'

    def _draw_thds(self, *events):
        if self.control_panels['ClassificationControl'].target in self.accepted_sources:
            spans = self.figure.select(tags='thd')
            spans.sort(key=lambda x: x.id)
            for i, span in enumerate(spans):
                if i < len(self.control_panels['ClassificationControl'].values):
                    span.location = self.control_panels['ClassificationControl'].values[i]
                    span.visible = self.control_panels['ClassificationControl'].show_thds
                else:
                    span.visible = False


class RateFigure(ThdLogFigure):
    title = 'Rates'
    accepted_sources = ['half-life', 'fit1', 'fit2', 'TF_rate']
    y_label = 'Rate (min⁻¹)'


class PFactFigure(ThdLogFigure):
    title = 'PFact'
    accepted_sources = ['pfact']  # list of names of sources which this plot accepts from parent controller
    y_label = 'Protection factor'


class FitResultFigure(BokehFigurePanel):
    title = 'Fit Result'
    accepted_sources = ['fr_pfact', 'uptake_corrected']
    y_label = 'Uptake corrected'
    x_label = 'Time'

    def __init__(self, parent, *args, **params):
        #todo refactor controllers to dict (Or better yet get them yourself from parent)
        super(FitResultFigure, self).__init__(parent, *args, **params)
        self.control_panels['FittingQuality'].param.watch(self._redraw_event, ['x_axis_type'])

    def _redraw_event(self, *events):
        self.redraw(x_axis_type=self.control_panels[0].x_axis_type.lower())

    def draw_figure(self, **kwargs):
        fig = super().draw_figure(x_axis_type=self.control_panels['FittingQuality'].x_axis_type.lower())
        return fig

    def render_sources(self, src_dict):
        for name, source in src_dict.items():
            func_name = 'line' if 'fr' in name else 'circle'  #todo default renderers
            glyph_func = getattr(self.figure, func_name)
            color = DEFAULT_COLORS[name]
            renderer = glyph_func(x='time', y='uptake', color=color, legend_label=name, source=source)
            self.renderers[name] = renderer

        if self.renderers:
            self.figure.legend.location = "bottom_right"


class ProteinFigure(FigurePanel):  #todo maybe it shouldnt be a figurepanel (it shoulnntr)
    title = 'Protein View'
    accepted_sources = ['pfact']

    js_files = {'ngl': "https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.37/dist/ngl.js"}

    def __init__(self, *args, **params):
        super(ProteinFigure, self).__init__(*args, **params)
        self.ngl_view = NGLViewer(sizing_mode='stretch_both')

        params = ['rcsb_id', 'no_coverage', 'representation', 'spin']
        self.parent.control_panels['ProteinViewControl'].param.watch(self._update_event, params)

    def _data_updated_callback(self, attr, old, new):
        self._update_colors(new['r_number'], new['color'])

    def render_sources(self, src_dict):
        for name, source in src_dict.items():
            self._update_colors(source.data['r_number'], source.data['color'])

    @property
    def no_coverage(self):
        return self.parent.control_panels['ProteinViewControl'].no_coverage

    def _update_colors(self, r_number, color_arr):
        r_start = r_number[0]
        color_list = list(color_arr)
        if r_start < 1:
            remove_num = 1 - r_start
            color_list = color_list[remove_num:]
        else:
            fill_num = r_start - 1
            color_list = fill_num*[self.parent.control_panels['ProteinViewControl'].no_coverage] + color_list

        self.ngl_view.color_list = color_list

    def _update_event(self, event):
        setattr(self.ngl_view, event.name, event.new)

    @property
    def panel(self):
        return self.ngl_view


class LogFigure(FigurePanel):
    title = 'Log'

    def __init__(self, *args, **params):
        super(LogFigure, self).__init__(*args, **params)
        self.markdown = LoggingMarkdown('### Log Window \n', sizing_mode='stretch_both')
        #todo config log level in options

        sh = logging.StreamHandler(self.markdown)
        sh.terminator = '  \n'
        formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', "%Y-%m-%d %H:%M:%S")
        sh.setFormatter(formatter)
        sh.setLevel(logging.DEBUG)
        self.parent.logger.addHandler(sh)

        #  todo add function on base class for doing these things (with try/except)
        self.parent.control_panels['OptionsPanel'].param.watch(self._update_log_level, ['log_level'])
        self.parent.control_panels['OptionsPanel'].param.trigger('log_level')

        # temporary to test logging
        self.parent.control_panels['DeveloperPanel'].param.watch(self._dev_btn, ['test_btn'])

    def _update_log_level(self, event):
        print('set log level', event.new)
        self.parent.logger.setLevel(event.new)

    def _dev_btn(self, events):
        # Temporary to test logging
        self.parent.logger.debug('dev test message')

    @property
    def panel(self):
        return self.markdown