import logging
import warnings

import panel as pn
import param


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

    _type = "base"

    sources = param.Dict(
        {}, doc="Dictionary of source objects available for plotting", precedence=-1
    )
    transforms = param.Dict({}, doc="Dictionary of transforms")
    opts = param.Dict({}, doc="Dictionary of formatting options (opts)")
    views = param.Dict({}, doc="Dictionary of views")

    loggers = param.Dict({}, doc="Dictionary of loggers")

    def __init__(self, control_panels, **params):
        super(MainController, self).__init__(**params)
        # self.client = client if client else Client()
        # check for client adress in config

        self.control_panels = {
            ctrl.name: ctrl(self) for ctrl in control_panels
        }  # todo as param?

        self.template = None  # Panel template (remove?)

        self.update()  # todo check to see if this is really needed

    # from lumen.target.Target
    def _rerender(self, *events, invalidate_cache=False):
        self._update_views(invalidate_cache=invalidate_cache)

    # todo remove?
    def _update_views(self, invalidate_cache=True, update_views=True, events=[]):
        warnings.warn("update view is deprecated", DeprecationWarning)
        for view in self.views.values():
            view.update()

    @property
    def panel(self):
        warnings.warn("panel property is deprecated", DeprecationWarning)
        # todo remove?
        return self.template

    def update(self):
        for view in self.views.values():
            view.update()


class PyHDXController(MainController):
    """
    Main controller for PyHDX web application.

    """

    _type = "pyhdx"

    sample_name = param.String(doc="Name describing the selected protein(s) state")

    def __init__(self, *args, **kwargs):
        super(PyHDXController, self).__init__(*args, **kwargs)

    @property
    def logger(self):
        return self.loggers["pyhdx"]
