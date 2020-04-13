
from .log import setup_custom_logger
from .panels import FileInputPanel, CoveragePanel
from pyhdx.pyhdx import PeptideCSVFile, KineticsSeries

logger = setup_custom_logger('root')
logger.debug('main message')


import param
import panel as pn
from jinja2 import Environment, FileSystemLoader
import holoviews as hv
import os

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
    peptides = param.ClassSelector(PeptideCSVFile)  #class with all peptides to be considered
    series = param.ClassSelector(KineticsSeries)

    def __init__(self, template, panels, **params):
        super(Controller, self).__init__(**params)
        template = env.get_template('template.html')

        tmpl = pn.Template(template=template)
     #   tmpl.nb_template.globals['get_id'] = make_globally_unique_id

        self.fileinput = FileInputPanel(self)
        self.coverage = CoveragePanel(self)

        # tmpl = pn.Template(template)
        tmpl.add_panel('controller', self.fileinput.control_panel)
        tmpl.add_panel('coverage', self.coverage.control_panel)
        tmpl.add_panel('scene3d', self.coverage.view_panel)
        tmpl.add_panel('slice_j', hv.Curve([1, 2, 3]))
        tmpl.add_panel('slice_k', hv.Curve([1, 2, 3]))


        #tmpl.add_panel('B', hv.Curve([1, 2, 3]))

        self.template = tmpl
      #  self.panels = [panel(self) for panel in panels]


    @property
    def servable(self):
        return self.template.servable

    @param.depends('data')
    def _test(self):
        print("hoi, data changed")


class HDXParams(param.Parameterized):
    """class for HDXPanel shared resources"""
    data = param.Array()

