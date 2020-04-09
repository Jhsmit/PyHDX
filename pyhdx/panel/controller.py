import param
import panel as pn
from jinja2 import Environment, FileSystemLoader
import holoviews as hv
import os

from bokeh.util.serialization import make_globally_unique_id


#from .panels import FileInputPanel
pth = os.path.dirname(__file__)

env = Environment(loader=FileSystemLoader(pth))

class Controller(param.Parameterized):
    """
    controller for main panels layout
    and has panels for each tabin the main layout

    """

    data = param.Array()


    def __init__(self, template, panels, **params):
        super(Controller, self).__init__(**params)
        template = env.get_template('template.html')

        tmpl = pn.Template(template=template)
     #   tmpl.nb_template.globals['get_id'] = make_globally_unique_id

        # tmpl = pn.Template(template)
        tmpl.add_panel('controller', hv.Curve([1, 2, 3]))
        tmpl.add_panel('scene3d', hv.Curve([1, 2, 3]))
        tmpl.add_panel('slice_i', hv.Curve([1, 2, 3]))
        tmpl.add_panel('slice_j', hv.Curve([1, 2, 3]))
        tmpl.add_panel('slice_k', hv.Curve([1, 2, 3]))


        #tmpl.add_panel('B', hv.Curve([1, 2, 3]))

        self.template = tmpl
      #  self.panels = [panel(self) for panel in panels]


    @property
    def servable(self):
        return self.template.servable

    @param.depends('hdxparams.data')
    def _test(self):
        print("hoi")


class HDXParams(param.Parameterized):
    """class for HDXPanel shared resources"""
    data = param.Array()

