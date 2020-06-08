from .base import FigurePanel, DEFAULT_RENDERERS, DEFAULT_COLORS
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column
from bokeh.models import LabelSet, ColumnDataSource, HoverTool, GlyphRenderer, Span
from bokeh.models.markers import Triangle, Circle, Diamond
import panel as pn
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

NGL_HTML = """
<div id="viewport" style="width:100%; height:100%;"></div>
<script>
stage = new NGL.Stage("viewport");
stage.loadFile("rcsb://1NKT.mmtf", {defaultRepresentation: true});
</script>
"""


class CoverageFigure(FigurePanel):

    def __init__(self, *args, **params):
        super(CoverageFigure, self).__init__(*args, **params)

        self.figures = [figure()]  #todo deprecate multiple figures
        self.figures[0].xaxis.axis_label = 'Residue number'

        self.layout = column(*self.figures, sizing_mode='stretch_both')
        self.label_set = LabelSet()
        self.bk_pane = pn.pane.Bokeh(self.layout, sizing_mode='stretch_both')

        #hook up main controller watchers
       # self.parent.param.watch(self._update, ['series'])

        self.ctrl = self.controllers[0]  # Only one controller for this Figure
        self.fit_ctrl = self.controllers[1]

        #hook up side watchers
        self.ctrl.param.watch(self._update_labels, ['labels'])
        self.ctrl.param.watch(self._update_index, ['index'])
        self.ctrl.param.watch(self._update, ['wrap', 'aa_per_subplot'])

        self.fit_ctrl.param.watch(self._draw_blocks, ['show_blocks', 'block_mode', 'max_combine', 'max_join',
                                                     'block_size', 'initial_block'])

    def _render_figure(self):
        raise NotImplementedError()

    @property
    def peptide_measurement(self):
        """return the peptide measurement currently plotted"""
        return self.parent.series[self.ctrl.index]

    @property
    def figure(self):
        return self.figures[0]

    def _draw_blocks(self, *events):
        for event in events:
            spans = self.figure.select(tags=['blocks'])
            if event.name == 'show_blocks' and spans:
               # if spans: # spans are already drawn, just toggle visible
                for span in spans:  # lovely span
                    span.visible = event.new
            else:  # All spans need to be (re)drawn
                for span in spans:
                    self.figure.center.remove(span)
                for pos in self.fit_ctrl.fit_block_edges:
                 #   print(pos)
                    span = Span(location=pos, dimension='height')
                    span.tags = ['blocks']
                    self.figure.add_layout(span)
                  #  print(self.figure.select(tags=['blocks']))

            #self.bk_pane.object = self.figure
            self.bk_pane.param.trigger('object')

    def _update(self, *events):
        print('series change on coverage fig')
        self.layout, self.figures, self.label_set = \
            _bokeh_coverage(self.peptide_measurement, self.ctrl.wrap, self.ctrl.aa_per_subplot)
        self.figures[-1].xaxis.axis_label = 'Residue number'

        self.label_set.visible = self.ctrl.labels
        self.bk_pane.object = self.layout

    def _update_labels(self, *events): #todo direct link?
        for event in events:
            self.label_set.visible = event.new

        self.bk_pane.param.trigger('object')

    def _update_index(self, *events):
        new_dict = {}
        new_dict['c'] = self._get_color()
        new_dict['uptake'] = self.peptide_measurement.data['uptake']
        new_dict['scores'] = self.peptide_measurement.data['scores']
        for fig in self.figures:
            fig.renderers[0].data_source.data.update(new_dict)
        self.bk_pane.param.trigger('object')

    def _get_color(self):
        #todo move function to pyhdx.plot
        #todo also update pallette in that function

        cmap = mpl.cm.get_cmap('jet')
        c_rgba = cmap(self.peptide_measurement.data['scores'] / 100)
        c = [mpl.colors.to_hex(color) for color in c_rgba]

        return list(c)


class RateFigure(FigurePanel):
    def __init__(self, *args, **params):
        super(RateFigure, self).__init__(*args, **params)

        self.figure = self.draw_figure()
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both')

        #todo refactor as kwargs?
        self.ctrl = self.controllers[1]  # classification controller
        self.fit_ctrl = self.controllers[0]
        self.ctrl.param.watch(self._draw_thds, ['values', 'show_thds'])
        self.parent.param.watch(self._update_rates, ['fit_results'])
        self.parent.param.watch(self._update_colors, ['rate_colors'])
        self.fit_ctrl.param.watch(self._draw_blocks, ['show_blocks', 'block_mode', 'max_combine', 'max_join',
                                                     'block_size', 'initial_block'])

    def draw_figure(self):
        """makes bokeh figure and returns it"""
        fig = figure(y_axis_type="log", tools='pan,wheel_zoom,box_zoom,save,reset,hover')
        fig.xaxis.axis_label = 'Residue number'
        fig.yaxis.axis_label = 'Rate'  # oh boy  #todo units?

        doc = curdoc()
        print(doc)

        #        DEFAULT_RENDERERS
        for k, v in DEFAULT_RENDERERS.items():
            glyph_func = getattr(fig, v)
            source = ColumnDataSource({name: [] for name in ['r_number', 'rate', 'color']})
            renderer = glyph_func(x='r_number', y='rate', color='color', source=source, legend_label=k, size=10,
                                  name=k)
            renderer.tags = ['rate']

        # spans for threshold lines
        for _ in range(self.controllers[1].param['num_classes'].bounds[1] - 1):  # todo refactor controller access
            sp = Span(location=0, dimension='width')
            sp.tags = ['thd']
            sp.visible = False
            fig.add_layout(sp)

        hover = fig.select(dict(type=HoverTool))
        hover.tooltips = [('Residue', '@r_number{int}'), ('Rate', '@rate')]
        hover.mode = 'vline'
        fig.legend.click_policy = 'hide'
        return fig

    def _update_rates(self, event):
        print('rates array update, renew')

        #todo maybe not if the user has already set it
        self.r_max = np.log(1 - 0.98) / - self.parent.series.times[1]

        #todo only redraw whats nessecary dependent on events?
        renderers = self.figure.select(tags='rate', type=GlyphRenderer)
        for renderer in renderers:
            array = self.parent.fit_results[renderer.name]['rates']
            new_dict = {name: array[name] for name in renderer.data_source.column_names if name in array.dtype.names}
            renderer.data_source.data.update(new_dict)

        self.bk_pane.param.trigger('object')

    def _update_colors(self, event):
        print('colorst in ratefigure')
        # #todo jslink colors to parent.rate_colors?
        # for event in events:
        renderers = self.figure.select(tags='rate', type=GlyphRenderer)
        for renderer in renderers:
            renderer.data_source.data.update({'color': self.parent.rate_colors[renderer.name]})
        self.bk_pane.param.trigger('object')

    def _draw_thds(self, *events):
        #todo check events and draw according to those? (events are triggers)
        spans = self.figure.select(tags='thd')
        spans.sort(key=lambda x: x.id)
        for i, span in enumerate(spans): # spanspanspam
            if i < len(self.ctrl.values):
                span.location = self.ctrl.values[i]
                span.visible = self.ctrl.show_thds
            else:
                span.visible = False


        # #remove everything
        # for some reason this wasnt a satisfactory solution
        # for span in self.figure.select(tags='thd'):
        #     self.figure.center.remove(span)
        # print("Values'", self.ctrl.values)
        # for value in self.ctrl.values:
        #     print('value in loop', value)
        #     span = Span(location=value, dimension='width')
        #     span.tags = ['thd']
        #     self.figure.add_layout(span)
        #     print(self.figure.center)

        self.bk_pane.param.trigger('object')

        #https://docs.bokeh.org/en/latest/docs/user_guide/layout.html
        #http://docs.bokeh.org/en/latest/docs/user_guide/annotations.html#spans
        #self.figure.add_layout()

    def _draw_blocks(self, *events):
        #duplicate code on CoverageFigure (also for ctrl reference)
        for event in events:
            spans = self.figure.select(tags=['blocks'])
            if event.name == 'show_blocks' and spans:
               # if spans: # spans are already drawn, just toggle visible
                for span in spans:  # lovely span
                    span.visible = event.new
            else:  # All spans need to be (re)drawn
                for span in spans:
                    self.figure.center.remove(span)
                for pos in self.fit_ctrl.fit_block_edges:
                 #   print(pos)
                    span = Span(location=pos, dimension='height')
                    span.tags = ['blocks']
                    self.figure.add_layout(span)
                  #  print(self.figure.select(tags=['blocks']))

            #self.bk_pane.object = self.figure
            self.bk_pane.param.trigger('object')

class FitResultFigure(FigurePanel):
    def __init__(self, *args, **params):
        super(FitResultFigure, self).__init__(*args, **params)

        self.ctrl = self.controllers[0]  #update to dict
        self.figure = self.draw_figure()
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both') #todo move to base class as property

        self.ctrl.param.watch(self._update_index, ['peptide_index'])
        self.ctrl.param.watch(self._redraw, ['x_axis_type'])
        self.parent.param.watch(self._update_fits, ['fit_results'])
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
        timepoints = self.parent.series.times
        datapoints = self.parent.series.scores_peptides.T[self.ctrl.peptide_index]
        #datapoints = [pm.data[self.ctrl.peptide_index]['scores'] for pm in self.parent.series]
        renderer = self.figure.select('Exp')[0]
        renderer.data_source.data.update({'time': timepoints, 'uptake': datapoints})
        self.bk_pane.param.trigger('object')

    @property
    def fit_timepoints(self):
        timepoints = self.parent.series.times
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


class ProteinFigure(FigurePanel):  #todo maybe it shouldnt be a figurepanel

    def __init__(self, *args, **params):
        super(ProteinFigure, self).__init__(*args, **params)

        self.html_panel = pn.pane.HTML(NGL_HTML, height=500, width=500)

    @property
    def panel(self):
        return pn.panel(self.html_panel)
