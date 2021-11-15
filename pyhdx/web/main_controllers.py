import logging
import param
import panel as pn

from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx import VERSION_STRING
from pyhdx.models import HDXMeasurementSet
from panel.template import BaseTemplate

from functools import partial
from dask.distributed import Client

from pyhdx.web.sources import PyHDXSource, AppSourceBase


class MainController(param.Parameterized):
    """
    Base class for application main controller
    Subclass to extend

    Parameters
    ----------
    control_panels : :obj:`list`
        List of strings referring to which ControlPanels to use for this MainController instance
        Should refer to subclasses of :class:`~pyhdx.panel.base.ControlPanel`
    client : dask client


    Attributes
    ----------

    doc : :class:`~bokeh.document.Document`
        Currently active Bokeh document
    logger : :class:`~logging.Logger`
        Logger instance
    control_panels : :obj:`dict`
        Dictionary with :class:`~pyhdx.panel.base.ControlPanel` instances (__name__ as keys)
    figure_panels : :obj:`dict`
        Dictionary with :class:`~pyhdx.panel.base.FigurePanel` instances (__name__ as keys)

    """

    _type = 'base'

    sources = param.Dict({}, doc='Dictionary of source objects available for plotting', precedence=-1)
    filters = param.Dict({}, doc="Dictionary of filters")
    opts = param.Dict({}, doc="Dictionary of formatting options (opts)")
    views = param.Dict({}, doc="Dictionary of views")

    logger = param.ClassSelector(logging.Logger, doc="Logger object")

    def __init__(self, control_panels, client=False, **params):
        super(MainController, self).__init__(**params)
        self.client = client if client else Client()
        if self.logger is None:
            self.logger = logging.getLogger(str(id(self)))

        self.control_panels = {ctrl.name: ctrl(self) for ctrl in control_panels}  #todo as param?

        self.template = None   # Panel template
        self.future_queue = []  # queue of tuples: (future, callback)

        self._update_views()
        self.start()

    # from lumen.target.Target
    def _rerender(self, *events, invalidate_cache=False):
        self._update_views(invalidate_cache=invalidate_cache)

    def _update_views(self, invalidate_cache=True, update_views=True, events=[]):
        for view in self.views.values():
            view.update()

    @property
    def panel(self):
        return self.template

    def update(self):
        for view in self.views:
            view.update()

    def check_futures(self):
        if self.future_queue:
            for future, callback in self.future_queue[:]:
                if future.status == 'finished':
                    callback(future)
                    self.future_queue.remove((future, callback))

    def start(self):
        refresh_rate = 1000
        pn.state.add_periodic_callback(
            self.check_futures, refresh_rate
        )


class PyHDXController(MainController):
    """
    Main controller for PyHDX web application.

    """

    _type = 'pyhdx'

    sample_name = param.String(doc='Name describing the selected protein(s) state')

    def __init__(self, *args, **kwargs):
        super(PyHDXController, self).__init__(*args, **kwargs)
    #
    # @param.depends('data_objects', watch=True)
    # def _datasets_updated(self):
    #     if len(self.data_objects) == 0:
    #         self.sample_name = ''
    #     elif len(self.data_objects) == 1:
    #         self.sample_name = str(next(iter(self.data_objects.keys())))
    #     elif len(self.data_objects) < 5:
    #         self.sample_name = ', '.join(self.data_objects.keys())
    #
    # @param.depends('sample_name', watch=True)
    # def _update_name(self):
    #     self.template.header[0].title = VERSION_STRING + ': ' + self.sample_name

