from __future__ import annotations

import logging
import warnings
from io import StringIO
from datetime import datetime
from typing import Optional, Tuple, List, Dict, Type, TYPE_CHECKING

import param
import yaml
from omegaconf import OmegaConf
from pyhdx.config import cfg
from pyhdx.support import clean_types
import pyhdx

if TYPE_CHECKING:
    from pyhdx.web.base import ControlPanel


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

    def __init__(self, control_panels: List[Tuple[Type[ControlPanel], Dict]], **params):
        super(MainController, self).__init__(**params)

        self.control_panels = {
            ctrl.name: ctrl(self, **kwargs) for ctrl, kwargs in control_panels
        }  # todo as param?

        self.template = None  # Panel template (remove?)

        self.session_time = datetime.now()
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

        self.log_io = StringIO()
        sh = logging.StreamHandler(self.log_io)
        sh.terminator = "  \n"
        # sh.setLevel(logging.CRITICAL)
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s]: %(message)s", "%Y-%m-%d %H:%M:%S"
        )
        sh.setFormatter(formatter)
        self.logger.addHandler(sh)

    @property
    def logger(self):
        return self.loggers["pyhdx"]

    def _get_file_header(self) -> str:
        """
        Returns header for txt file outputs with version and date/time information

        """

        s = f"# {pyhdx.VERSION_STRING} \n"
        s += (
            f"# {self.session_time.strftime('%Y/%m/%d %H:%M:%S')}"
            f"({int(self.session_time.timestamp())})\n"
        )

        return s

    def hdx_spec_callback(self) -> Optional[StringIO]:
        """
        Get a StringIO with input HDX measurement specifications.

        """
        input_controllers = {"PeptideFileInputControl", "PeptideRFUFileInputControl"}
        input_ctrls = self.control_panels.keys() & input_controllers
        if len(input_ctrls) == 1:
            input_ctrl = self.control_panels[list(input_ctrls)[0]]

            s = yaml.dump(clean_types(input_ctrl.hdx_spec), sort_keys=False)
            output = self._get_file_header() + "\n" + s
            sio = StringIO(output)

            return sio
        else:
            return None

    def config_callback(self) -> StringIO:
        """
        Get a StringIO with global configuration settings.

        """

        masked_conf = OmegaConf.masked_copy(cfg.conf, cfg.conf.keys() - {"server"})
        s = OmegaConf.to_yaml(masked_conf)

        output = self._get_file_header() + "\n" + s
        sio = StringIO(output)

        return sio

    def user_settings_callback(self) -> StringIO:
        """
        Get a StringIO with user settings.

        """
        user_dict = self.sources["metadata"].get("user_settings")
        s = yaml.dump(clean_types(user_dict), sort_keys=False)

        output = self._get_file_header() + "\n" + s
        sio = StringIO(output)

        return sio

    def log_callback(self) -> StringIO:
        """
        Get a StringIO with the full log.

        """

        self.log_io.seek(0)
        s = self.log_io.read()

        output = self._get_file_header() + "\n" + s
        sio = StringIO(output)

        return sio


# single amide slider only first?
class PeptideController(MainController):

    """Object which models D-uptake of the peptide in time"""

    _type = "peptide"

    def __init__(self, *args, **params) -> None:
        super().__init__(*args, **params)

        # source = TableSource()
        # self.sources['main'] = source
        # view = hvCurveView(source=source)
        # self.views['peptide'] = view

    # def _action_reload(self):
    #     self.model = PeptideUptakeModel(self.fasta_sequence, self.temperature, self.pH)
    #
    # def get_plot(self):
    #     func = partial(hv.Curve, kdims='index', vdims='p_ND')
    #
    #     dm = hv.DynamicMap(func, streams=[stream])
    #     dm = dm.apply.opts(logx=True)
