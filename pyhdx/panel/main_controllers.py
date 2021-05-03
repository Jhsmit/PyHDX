import logging
import param
import panel as pn

from pyhdx.models import PeptideMasterTable, KineticsSeries
from pyhdx import VERSION_STRING_SHORT

from panel.template import BaseTemplate


from lumen.sources import Source
from lumen.filters import FacetFilter

from functools import partial

class MainController(param.Parameterized):
    """
    Base class for application main controller
    Subclass to extend

    Parameters
    ----------
    control_panels : :obj:`list`
        List of strings referring to which ControlPanels to use for this MainController instance
        Should refer to subclasses of :class:`~pyhdx.panel.base.ControlPanel`
    cluster : :obj:`str`
        IP:port address for Dask cluster (optional)

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
    sources = param.Dict({}, doc='Dictionary of source objects available for plotting', precedence=-1)
    transforms = param.Dict({}, doc='Dictionary of transforms')
    filters = param.Dict({}, doc="Dictionary of filters")
    opts = param.Dict({}, doc="Dictionary of formatting options (opts)")
    views = param.Dict({}, doc="Dictionary of views")

    logger = param.ClassSelector(logging.Logger, doc="Logger object")

    def __init__(self, control_panels, cluster=None, **params):
        super(MainController, self).__init__(**params)
        self.cluster = cluster
        self._doc = pn.state.curdoc
        if self.logger is None:
            self.logger = logging.getLogger(str(id(self)))

        self.control_panels = {ctrl.name: ctrl(self) for ctrl in control_panels}  #todo as param?

        self.template = None   # Panel template

        for filt in self.filters.values():
            if isinstance(filt, FacetFilter):
                continue
            print(filt)
            filt.param.watch(partial(self._rerender, invalidate_cache=True), 'value')

        for trs in self.transforms.values():
            if hasattr(trs, 'updated'):
                trs.param.watch(partial(self._rerender, invalidate_cache=True), 'updated')

        self._update_views()

    # from lumen.target.Target
    def _rerender(self, *events, invalidate_cache=False):
        self._update_views(invalidate_cache=invalidate_cache)

    def _update_views(self, invalidate_cache=True, update_views=True, events=[]):
        for view in self.views.values():
            view.update(invalidate_cache=invalidate_cache)

    @property
    def doc(self):
        """ :class:`~bokeh.document.document.Document`: Bokeh document for the application"""
        return self._doc or pn.state.curdoc

    def __panel__(self):
        # This does something but not as expected
        return self.template

    @property
    def panel(self):
        return self.template

    def update(self):
        for view in self.views:
            view.update()

    def start(self):
        refresh_rate = 50
        self._cb = pn.state.add_periodic_callback(
            self.update, refresh_rate
        )


class PyHDXController(MainController):
    """
    Main controller for PyHDX web application.

    """

    fit_objects = param.Dict(default={}, doc='Dictionary for all datasets (KineticsFitting objects)')

    # for guesses (nested): <fit name>: {state1: state2:, ...
    # for global fit (batch): <fit name>: fit_result_object
    # for global fit (series): <fit name>: {state1: <fit_result_object>, state2:....}
    fit_results = param.Dict({}, doc='Dictionary of fit results', precedence=-1)

    sample_name = param.String(doc='Name describing the selected protein state')

    def __init__(self, *args, **kwargs):
        super(PyHDXController, self).__init__(*args, **kwargs)

    @param.depends('fit_objects', watch=True)
    def _datasets_updated(self):
        if len(self.fit_objects) == 0:
            self.sample_name = ''
        elif len(self.fit_objects) == 1:
            self.sample_name = str(next(iter(self.fit_objects.keys())))
        elif len(self.fit_objects) < 5:
            self.sample_name = ', '.join(self.fit_objects.keys())

    @param.depends('sample_name', watch=True)
    def _update_name(self):
        self.template.header[0].title = VERSION_STRING_SHORT + ': ' + self.sample_name


class ComparisonController(MainController):
    """
    Main controller for binary comparison web application.
    """

    datasets = param.Dict(default={}, doc='Dictionary for all datasets')
    comparisons = param.Dict(default={}, doc='Dictionary for all comparisons (should be in sources)')
