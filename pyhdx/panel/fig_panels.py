from .base import FigurePanel, DEFAULT_RENDERERS
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure
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

        self.figures = [figure()]
        self.figures[0].xaxis.axis_label = 'Residue number'

        self.layout = column(*self.figures, sizing_mode='stretch_both')
        self.label_set = LabelSet()
        self.bk_pane = pn.pane.Bokeh(self.layout, sizing_mode='stretch_both')

        #hook up main controller watchers
       # self.parent.param.watch(self._update, ['series'])

        self.ctrl = self.controllers[0]  # Only one controller for this Figure

        #hook up side watchers
        self.ctrl.param.watch(self._update_labels, ['labels'])
        self.ctrl.param.watch(self._update_index, ['index'])
        self.ctrl.param.watch(self._update, ['wrap', 'aa_per_subplot'])

    def _render_figure(self):
        raise NotImplementedError()

    @property
    def peptide_measurement(self):
        """return the peptide measurement currently plotted"""
        return self.parent.series[self.ctrl.index]

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

        self.figure = figure(y_axis_type="log", tools='pan,wheel_zoom,box_zoom,save,reset,hover')
        self.figure.xaxis.axis_label = 'Residue number'
        self.figure.yaxis.axis_label = 'Rate (min⁻¹)'  # oh boy

#        DEFAULT_RENDERERS
        for k, v in DEFAULT_RENDERERS.items():
            glyph_func = getattr(self.figure, v)
            source = ColumnDataSource({name: [] for name in ['r_number', 'rate', 'color']})
            renderer = glyph_func(x='r_number', y='rate', color='color', source=source, legend_label=k, size=10,
                                  name=k)
            renderer.tags = ['rate']

        #spans for threshold lines
        for _ in range(self.controllers[1].param['num_classes'].bounds[1] - 1):  # todo refactor controller access
            sp = Span(location=0, dimension='width')
            sp.tags = ['thd']
            sp.visible = False
            self.figure.add_layout(sp)

        hover = self.figure.select(dict(type=HoverTool))
        hover.tooltips = [('Residue', '@r_number{int}'), ('Rate', '@rate')]
        hover.mode = 'vline'
        self.figure.legend.click_policy = 'hide'
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both')

        #todo refactor as kwargs?
        self.ctrl = self.controllers[1]  # classification controller
        self.ctrl.param.watch(self._draw_thds, ['values', 'show_thds'])
        self.parent.param.watch(self._update_rates, ['fit_results'])
        self.parent.param.watch(self._update_colors, ['rate_colors'])

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


class ProteinFigure(FigurePanel):  #todo maybe it shouldnt be a figurepanel

    def __init__(self, *args, **params):
        super(ProteinFigure, self).__init__(*args, **params)

        self.html_panel = pn.pane.HTML(NGL_HTML, height=500, width=500)

    @property
    def panel(self):
        return pn.panel(self.html_panel)
