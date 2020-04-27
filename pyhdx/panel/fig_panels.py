from .base import FigurePanel
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import LabelSet, ColumnDataSource, HoverTool
from bokeh.models.markers import Triangle, Circle, Diamond
import panel as pn
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt


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

        hover = self.figure.select(dict(type=HoverTool))
        hover.tooltips = [('Residue', '@r_number{int}'), ('Rate', '@rate')]
        hover.mode = 'vline'
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both')

        self.fit_renderers = {}
        self.line_renderers = {}

        #todo refactor as kwargs?
        self.ctrl = self.controllers[1]  # classification controller
        self.ctrl.param.watch(self._draw_thds, ['values'])
        self.parent.param.watch(self._renew, ['fit_results'])
        self.parent.param.watch(self._update, ['series'])

    def _update(self, *events):
        #redraw plot because of new series

        DEFAULT_RENDERERS = {'fit1': 'triangle', 'fit2': 'circle'}
        DEFAULT_COLORS = {'fit1': 'blue', 'fit2': 'red'}

        self.fit_renderers = {}
        for k, v in self.parent.fit_results.items():
            array = v['rates']
            glyph_func = getattr(self.figure, DEFAULT_RENDERERS.get(k, 'diamond'))
            color = DEFAULT_COLORS.get(k, 'black')
            print([color]*len(array))
            source = ColumnDataSource({'r_number': array['r_number'], 'rate': array['rate'], 'color': [color]})  #todo add color
            renderer = glyph_func(x='r_number', y='rate', source=source, color='color', legend_label='asdf', size=10)

            #r = self.figure.triangle(x='r_number', y='rate', legend_label=k, source=source)#, color='fit1_color')
            #legend_label=k,
            #renderer = self.figure.add_glyph(source,  glyph=glyph)
            self.fit_renderers[k] = renderer

        self.figure.legend.click_policy = 'hide'
        #self.figure.circle(x='r_number', y='fit1_r1', legend_label='Fit 1 r1', source=source, color='green')
        #self.figure.circle(x='r_number', y='fit1_r2', legend_label='Fit 1 r2', source=source, color='yellow')

        self.bk_pane.param.trigger('object')

    def _renew(self, event):
        print('rates array update, renew')

        #todo maybe not if the user has already set it
        self.r_max = np.log(1 - 0.98) / - self.parent.series.times[1]

        #todo only redraw whats nessecary dependent on events?
        #new_dict = {name: self.parent.rates[name] for name in self.parent.rates.dtype.names}

        #todo this should be down with plot.select and naming the glyps
        # >> > plot.circle([1, 2, 3], [4, 5, 6], name="temp")
        # >> > plot.select(name="temp")
        # [GlyphRenderer(id='399d53f5-73e9-44d9-9527-544b761c7705', ...)]

        # >> > r = plot.circle([1, 2, 3], [4, 5, 6])
        # >> > r.tags = ["foo", 10]
        # >> > plot.select(tags=['foo', 10])
        # [GlyphRenderer(id='1de4c3df-a83d-480a-899b-fb263d3d5dd9', ...)]

        for key, renderer in self.fit_renderers.items():
            array = self.parent.fit_results[key]['rates']
            new_dict = {name: array[name] for name in renderer.data_source.column_names if name in array.dtype.names}

            #Temporary hack
            new_dict['color'] = [renderer.data_source.data['color'][0]]*len(array)
            renderer.data_source.data.update(new_dict)

        self.bk_pane.param.trigger('object')

    def _draw_thds(self, *events):
        #todo check events and draw according to those?
        print('draw thresholds')
        for event in events:
            print(event)

        #https://docs.bokeh.org/en/latest/docs/user_guide/layout.html
        #http://docs.bokeh.org/en/latest/docs/user_guide/annotations.html#spans
        #self.figure.add_layout()
