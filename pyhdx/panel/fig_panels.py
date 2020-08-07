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

NGL_HTML = """
<div id="viewport" style="width:100%; height:100%;"></div>
<script>
stage = new NGL.Stage("viewport");
stage.loadFile("rcsb://1NKT.mmtf", {defaultRepresentation: true});
</script>
"""


class CoverageFigure(FigurePanel):
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
            # with pn.io.unlocked():
            #     self.figure.x_range.start = source.data['start'].min() - 50
            #     self.figure.x_range.end = source.data['end'].max() + 50
            #     self.update()


class ThdLogFigure(FigurePanel):
    """base class for pfact / rates figure panels which both feature y log axis and thresholding"""

    def __init__(self, parent, controllers, *args, **params):
        #todo refactor controllers to dict (Or better yet get them yourself from parent)
        self.ctrl = controllers[1]  # classification controller
        self.fit_ctrl = controllers[0]
        super(ThdLogFigure, self).__init__(parent, controllers, *args, **params)

        self.ctrl.param.watch(self._draw_thds, ['values', 'show_thds'])

    def draw_figure(self):
        fig = figure(y_axis_type='log', tools='pan,wheel_zoom,box_zoom,save,reset')
        fig.min_border_left = MIN_BORDER_LEFT
        fig.xaxis.axis_label = 'Residue number'
        fig.yaxis.axis_label = self.y_label

        for _ in range(self.controllers[1].param['num_classes'].bounds[1] - 1):  # todo refactor controller access
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
    accepted_sources = ['half-life', 'fit1', 'fit2', 'TF_rate']
    y_label = 'Rate (min⁻¹)'


class PFactFigure(ThdLogFigure):
    accepted_sources = ['pfact']  # list of names of sources which this plot accepts from parent controller
    y_label = 'Protection factor'


class FitResultFigure(FigurePanel):
    accepted_sources = ['fr_pfact', 'uptake_corrected']
    y_label = 'Uptake corrected'
    x_label = 'Time'

    def __init__(self, parent, controllers, *args, **params):
        #todo refactor controllers to dict (Or better yet get them yourself from parent)
        super(FitResultFigure, self).__init__(parent, controllers, *args, **params)

        self.controllers[0].param.watch(self._redraw_event, ['x_axis_type'])

    def _redraw_event(self, *events):
        self.redraw(x_axis_type=self.controllers[0].x_axis_type.lower())

    def draw_figure(self, **kwargs):
        fig = super().draw_figure(x_axis_type=self.controllers[0].x_axis_type.lower())

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

class FitResultFigureOld(FigurePanelOld):
    def __init__(self, *args, **params):
        super(FitResultFigureOld, self).__init__(*args, **params)

        self.ctrl = self.controllers[0]  #update to dict
        self.figure = self.draw_figure()
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both') #todo move to base class as property

        self.ctrl.param.watch(self._update_index, ['peptide_index'])
        self.ctrl.param.watch(self._redraw, ['x_axis_type'])
        #self.parent.param.watch(self._update_fits, ['fit_results'])
        self.parent.param.watch(self._update_peptides, ['series'])
        #self.parent.param.watch(self._update_colors, ['rate_colors'])

        self.fit_data = {} # Dictionary with uptake values for peptides from fits

    def _redraw(self, *events):
        print('redraw')
        print(self.ctrl.x_axis_type.lower())
        self.figure = self.draw_figure()

        self._update_fits()
        self._update_peptides()
        self.bk_pane.object = self.figure

    def draw_figure(self):
        fig = figure(x_axis_type=self.ctrl.x_axis_type.lower())
        fig.axis.axis_label = 'Time'
        doc = curdoc()
        print(doc)
        renderers = ['fit1', 'fit2']
        for r in renderers:
            source = ColumnDataSource({'time': [], 'uptake': []})
            fig.line(x='time', y='uptake', color=DEFAULT_COLORS[r], source=source, legend_label=r, name=r)

        source = ColumnDataSource({'time': [], 'uptake': []})
        fig.circle(x='time', y='uptake', color='black', source=source, legend_label='Exp', name='Exp')
        return fig

    def _update_peptides(self, *events):
        timepoints = self.parent.series.timepoints
        datapoints = self.parent.series.scores_peptides.T[self.ctrl.peptide_index]
        #datapoints = [pm.data[self.ctrl.peptide_index]['scores'] for pm in self.parent.series]
        renderer = self.figure.select('Exp')[0]
        renderer.data_source.data.update({'time': timepoints, 'uptake': datapoints})
        self.bk_pane.param.trigger('object')

    @property
    def fit_timepoints(self):
        timepoints = self.parent.series.timepoints
        if self.ctrl.x_axis_type == 'Linear':
            time = np.linspace(0, timepoints.max(), num=250)
        elif self.ctrl.x_axis_type == 'Log':
            time = np.logspace(-2, np.log10(timepoints.max()), num=250)
        return time

    def _update_fits(self, *events):
        for k, v in self.parent.fit_results.items():
            fit_result = v['fitresult']
            if fit_result is None:
                continue
            #Issue #35

            if k == 'fit1':
                d_list = []
                for time in self.fit_timepoints:
                    p = fit_result.get_p(time)
                    p = np.nan_to_num(p)
                    d = self.parent.series.cov.X.dot(p)
                    d_list.append(d)
                d = np.vstack(d_list)
            elif k == 'fit2':
                d_list = []
                for time in self.fit_timepoints:
                    d = fit_result.get_d(time)
                    d_list.append(d)
                d = np.vstack(d_list)

            self.fit_data[k] = d
        self._update_index()

    def _update_index(self, *events):
        for k, v in self.fit_data.items():
            renderer = self.figure.select(k)[0]
            renderer.data_source.data.update({'time': self.fit_timepoints, 'uptake': v[:, self.ctrl.peptide_index]})

        #update measured points
        print('series', self.parent.series)
        datapoints = [pm.data[self.ctrl.peptide_index]['scores'] for pm in self.parent.series]
        renderer = self.figure.select('Exp')[0]
        renderer.data_source.data.update({'uptake': datapoints})

        self.bk_pane.param.trigger('object')


class ProteinFigure(FigurePanelOld):  #todo maybe it shouldnt be a figurepanel (it shoulnntr)

    def __init__(self, *args, **params):
        super(ProteinFigure, self).__init__(*args, **params)

        self.html_panel = pn.pane.HTML(NGL_HTML, height=500, width=500)

    @property
    def panel(self):
        #view = nglview.show_file(r"C:\Users\jhsmi\pp\pyHDX_paper\fig4_real_proteins\structures\hPREP_2.pdb")
        return None
        #return pn.panel(view)
        #return pn.panel(self.html_panel)
