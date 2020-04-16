
from .log import setup_custom_logger
from .base import ControlPanel
from .panels import FileInputPanel, CoveragePanel, RateConstantPanel
from .fig_panels import CoverageFigure
from pyhdx.pyhdx import PeptideCSVFile, KineticsSeries
from pyhdx.fitting import KineticsFitting

logger = setup_custom_logger('root')
logger.debug('main message')


import param
import panel as pn
from jinja2 import Environment, FileSystemLoader
import holoviews as hv
import os
import numpy as np

import matplotlib
matplotlib.use('agg') # for panel mpl support

from bokeh.util.serialization import make_globally_unique_id


pth = os.path.dirname(__file__)

env = Environment(loader=FileSystemLoader(pth))


class Controller(param.Parameterized):
    """
    controller for main panels layout
    and has panels for each tabin the main layout

    """

    data = param.Array()  # might not be needed, in favour of peptides
    rates = param.Array(doc='Output rates data')
    peptides = param.ClassSelector(PeptideCSVFile)  #class with all peptides to be considered
    series = param.ClassSelector(KineticsSeries)
    fitting = param.ClassSelector(KineticsFitting)

    def __init__(self, template, panels, **params):
        super(Controller, self).__init__(**params)
        template = env.get_template('template.html')

        tmpl = pn.Template(template=template)
     #   tmpl.nb_template.globals['get_id'] = make_globally_unique_id


        # Controllers
        self.fileinput = FileInputPanel(self)
        self.coverage = CoverageControl(self)#CoveragePanel(self)
        self.rate_panel = RateConstantPanel(self)


        #Figures
        self.coverage_figure = CoverageFigure(self, [self.coverage])  #parent, [controllers]

        # tmpl = pn.Template(template)
        tmpl.add_panel('input', self.fileinput.control_panel)
        tmpl.add_panel('coverage', self.coverage.panel)
        tmpl.add_panel('fitting', self.rate_panel.control_panel)
        tmpl.add_panel('scene3d', self.coverage_figure.panel)
        tmpl.add_panel('slice_j', self.rate_panel.view_panel)
        tmpl.add_panel('slice_k', hv.Curve([1, 2, 3]))


        #tmpl.add_panel('B', hv.Curve([1, 2, 3]))

        self.template = tmpl
      #  self.panels = [panel(self) for panel in panels]

    @param.depends('series', watch=True)
    def _series_changed(self):
        # This is triggered if the fileinput child panel yields a new KineticSeries
        print('series changed')

        self.fitting = KineticsFitting(self.series)
        #todo add errors here
        rate_fields = ['fit1', 'fit1_r1', 'fit1_r2', 'fit2', 'fit2_r1', 'fit2_r2']
        rates = np.zeros(self.series.cov.prot_len,
                              dtype=[('r_number', int)] + [(name, float ) for name in rate_fields])
        rates['r_number'] = self.series.cov.r_number

        self.rates = rates  # this assignement triggers downstream watchers

    @property
    def servable(self):
        return self.template.servable

    @param.depends('data')
    def _test(self):
        print("hoi, data changed")


class CoverageControl(ControlPanel):
    wrap = param.Integer(25, bounds=(0, None), doc='Number of peptides vertically before moving to the next row') # todo auto?
    aa_per_subplot = param.Integer(100, label='Amino acids per subplot')
    update = param.Action(lambda self: self.param.trigger('update'), label='Update')  #Triggers redraw of the figure (might not need this in favour of directly triggering from wrap and aa per subplot
    labels = param.Boolean(False, label='Labels')
    index = param.Integer(0, bounds=(0, None), doc='Current index of coverage plot in time')

    def __init__(self, *args, **params):
        super(CoverageControl, self).__init__(*args, **params)
        self.exposure_str = pn.widgets.StaticText(name='Exposure', value='0') # todo update to some param?
        self.parent.param.watch(self._update_series, ['series'])

    @property
    def peptide_measurement(self):
        return self.parent.series[self.index]

    def _update_series(self, event):
        print('coverage new series update index bounds')
        self.param['index'].bounds = (0, len(event.new) - 1)
        self.exposure_str.value = str(self.peptide_measurement.exposure)

    @param.depends('index', watch=True)
    def _update_index(self):
        self.exposure_str.value = str(self.peptide_measurement.exposure)

    @property
    def panel(self):
        col = pn.Column(self.param, self.exposure_str)
        #p = pn.Param(self.param, widgets={'file': pn.widgets.FileInput}) for exposure
        return pn.WidgetBox(col, pn.layout.VSpacer(), css_classes=['widget-box', 'custom-wbox'],
                            sizing_mode='stretch_height')
