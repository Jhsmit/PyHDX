from .base import FigurePanelOld, FigurePanel, DEFAULT_RENDERERS, DEFAULT_COLORS, MIN_BORDER_LEFT
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column
from bokeh.models import LabelSet, ColumnDataSource, HoverTool, GlyphRenderer, Span, Rect, Range1d
from bokeh.models.markers import Triangle, Circle, Diamond
import panel as pn
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import nglview

import param

NGL_HTML = """
<div id="viewport" style="width:100%; height:100%;"></div>
<script>
stage = new NGL.Stage("viewport");
stage.loadFile("rcsb://1NKT.mmtf", {defaultRepresentation: true});
</script>
"""


class CoverageFigure(FigurePanel):
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


class ThdLogFigure(FigurePanel):
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
        # todo check events and draw according to those? (events are triggers)
        if self.ctrl.target in self.accepted_sources:
            spans = self.figure.select(tags='thd')
            spans.sort(key=lambda x: x.id)
            for i, span in enumerate(spans): # spanspanspam
                if i < len(self.ctrl.values):
                    span.location = self.ctrl.values[i]
                    span.visible = self.ctrl.show_thds
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


class FitResultFigure(FigurePanel):
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


class ProteinFigure(FigurePanelOld):  #todo maybe it shouldnt be a figurepanel (it shoulnntr)

    def __init__(self, *args, **params):
        super(ProteinFigure, self).__init__(*args, **params)

        self.html_panel = pn.pane.HTML(NGL_HTML, height=500, width=500)

    @property
    def panel(self):
        view = nglview.show_file(r"C:\Users\jhsmi\pp\pyHDX_paper\v1\fig4_real_proteins\structures\hPREP.pdb")
        return pn.pane.IPyWidget(view)
        #return pn.panel(view)
        #return pn.panel(self.html_panel)
