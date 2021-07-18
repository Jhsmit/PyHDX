import logging
import param
import panel as pn

from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx import VERSION_STRING
from pyhdx.models import HDXMeasurementSet
from panel.template import BaseTemplate


from lumen.sources import Source
from lumen.filters import FacetFilter

from functools import partial
from dask.distributed import Client


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
    sources = param.Dict({}, doc='Dictionary of source objects available for plotting', precedence=-1)
    transforms = param.Dict({}, doc='Dictionary of transforms')
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

        for filt in self.filters.values():
            if isinstance(filt, FacetFilter):
                continue
            filt.param.watch(partial(self._rerender, invalidate_cache=True), 'value')

        for trs in self.transforms.values():
            if hasattr(trs, 'updated'):
                trs.param.watch(partial(self._rerender, invalidate_cache=True), 'updated')

        self._update_views()
        self.start()

    # from lumen.target.Target
    def _rerender(self, *events, invalidate_cache=False):
        self._update_views(invalidate_cache=invalidate_cache)

    def _update_views(self, invalidate_cache=True, update_views=True, events=[]):
        for view in self.views.values():
            view.update(invalidate_cache=invalidate_cache)

    def __panel__(self):
        # This does something but not as expected
        return self.template

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

    data_objects = param.Dict(default={}, doc='Dictionary for all datasets (HDXMeasurement objects)') # todo refactor

    # for guesses (nested): <fit name>: {state1: state2:, ...
    # for global fit (batch): <fit name>: fit_result_object
    # for global fit (series): <fit name>: {state1: <fit_result_object>, state2:....}
    fit_results = param.Dict({}, doc='Dictionary of fit results', precedence=-1)

    sample_name = param.String(doc='Name describing the selected protein(s) state')

    def __init__(self, *args, **kwargs):
        super(PyHDXController, self).__init__(*args, **kwargs)

    @param.depends('data_objects', watch=True)
    def _datasets_updated(self):
        if len(self.data_objects) == 0:
            self.sample_name = ''
        elif len(self.data_objects) == 1:
            self.sample_name = str(next(iter(self.data_objects.keys())))
        elif len(self.data_objects) < 5:
            self.sample_name = ', '.join(self.data_objects.keys())

    @param.depends('sample_name', watch=True)
    def _update_name(self):
        self.template.header[0].title = VERSION_STRING + ': ' + self.sample_name

    @property
    def hdx_set(self):
        """Returns combined HDXMeasurementSet of all currently added data objects"""
        #todo when alignments are added in, update this as (fixed) attribute

        return HDXMeasurementSet(list(self.data_objects.values()))


class ComparisonController(MainController):
    """
    Main controller for binary comparison web application.
    """

    datasets = param.Dict(default={}, doc='Dictionary for all datasets')
    comparisons = param.Dict(default={}, doc='Dictionary for all comparisons (should be in sources)')
