from .base import FigurePanel
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import LabelSet, ColumnDataSource
import panel as pn
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt


class CoverageFigure(FigurePanel):

    def __init__(self, *args, **params):
        super(CoverageFigure, self).__init__(*args, **params)

        self.figures = [figure()]
        self.layout = column(*self.figures, sizing_mode='stretch_both')
        self.label_set = LabelSet()
        self.bk_pane = pn.pane.Bokeh(self.layout, sizing_mode='stretch_both')

        #hook up main controller watchers
        self.parent.param.watch(self._update, ['series'])

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
        self.label_set.visible = self.ctrl.labels
        self.bk_pane.object = self.layout

    def _update_labels(self, *events): #todo direct link?
        for event in events:
            self.label_set.visible = event.new

        self.bk_pane.param.trigger('object')

    def _update_index(self, *events):
        color = self._get_color()
        for fig in self.figures:
            fig.renderers[0].data_source.data.update({'c': color})
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

        self.figure = figure(y_axis_type="log")
        self.figure.xaxis.axis_label = 'Residue number'
        self.figure.yaxis.axis_label = 'Rate (min⁻¹)'  # oh boy
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both')

        self.parent.param.watch(self._renew, ['rates'])
        self.parent.param.watch(self._update, ['series'])

        self.fit_renderers = []
        self.line_renderers = []

        #todo refactor as kwargs?
        self.ctrl = self.controllers[1]  # classification controller
        self.ctrl.param.watch(self._draw_thds, ['values'])

    def _update(self, *events):
        #draw plot because of new series

        source = ColumnDataSource({name: self.parent.rates[name] for name in self.parent.rates.dtype.names})

        self.line_renderers = []
        r = self.figure.triangle(x='r_number', y='fit1', legend_label='Fit 1', source=source, color='fit1_color')
        self.line_renderers.append(r)
        r = self.figure.circle(x='r_number', y='fit2', legend_label='Fit 2', source=source, color='fit2_color')
        self.line_renderers.append(r)
        self.figure.legend.click_policy = 'hide'

        #self.figure.circle(x='r_number', y='fit1_r1', legend_label='Fit 1 r1', source=source, color='green')
        #self.figure.circle(x='r_number', y='fit1_r2', legend_label='Fit 1 r2', source=source, color='yellow')

        self.bk_pane.param.trigger('object')

    def _renew(self, event):
        print('rates array update, renew')

        #todo maybe not if the user has already set it
        self.r_max = np.log(1 - 0.98) / - self.parent.series.times[1]

        new_dict = {name: self.parent.rates[name] for name in self.parent.rates.dtype.names}
        for renderer in self.fit_renderers:
            renderer.data_source.data.update(new_dict)

        self.bk_pane.param.trigger('object')

    def _draw_thds(self, *events):
        #todo check events and draw according to those?
        print('draw thresholds')

        #https://docs.bokeh.org/en/latest/docs/user_guide/layout.html
        #http://docs.bokeh.org/en/latest/docs/user_guide/annotations.html#spans
        #self.figure.add_layout()
