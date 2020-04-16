from .base import FigurePanel
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import LabelSet
import panel as pn

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
        self.ctrl.param.watch(self._update, ['update'])

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

    def _update_labels(self, event):
        self.label_set.visible = event.new
        self.bk_pane.param.trigger('object')

    def _update_index(self, event):
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
