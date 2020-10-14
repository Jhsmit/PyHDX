import logging
import param
import panel as pn

from pyhdx.models import PeptideMasterTable, KineticsSeries


class MainController(param.Parameterized):
    """
    Base class for application main controller
    Subclass to extend

    Parameters
    ----------
    control_panels: :obj:`list`
        List of strings referring to which ControlPanels to use for this MainController instance
        Should refer to subclasses of :class:`~pyhdx.panel.base.ControlPanel`
    figure_panels: :obj:`list`
        List of string referring to which FigurePanels to use for this MainController instance
        Should refer to subclasses :class:`~pyhdx.panel.base.FigurePanel`
    cluster: :obj:`str`
        IP:port address for Dask cluster (optional)

    Attributes
    ----------

    doc : :class:`~bokeh.document.Document`
        Currently active Bokeh document
    logger : :class:`~logging.Logger`
        Logger instance
    control_panels : :obj:`dict`
        Dictionary with :class:`~pyhdx.panel.base.ControlPanel` instances (__name__ as keys)
    figure_panels : :obj`dict`
        Dictionary with :class:`~pyhdx.panel.base.FigurePanel` instances (__name__ as keys)

    """
    sources = param.Dict({}, doc='Dictionary of ColumnDataSources available for plotting', precedence=-1)

    def __init__(self, control_panels, figure_panels, cluster=None, **params):
        super(MainController, self).__init__(**params)
        self.cluster = cluster
        self._doc = pn.state.curdoc
        self.logger = logging.getLogger(str(id(self)))

        #available_controllers = {cls.__name__: cls for cls in gen_subclasses(ControlPanel)}
        self.control_panels = {ctrl.name: ctrl(self) for ctrl in control_panels}

        #available_figures = {cls.__name__: cls for cls in gen_subclasses(FigurePanel)}
        self.figure_panels = {ctrl.name: ctrl(self) for ctrl in figure_panels}

    @property
    def doc(self):
        """ :class:`~bokeh.document.document.Document`: Bokeh document for the application"""
        return self._doc or pn.state.curdoc

    def publish_data(self, name, data_source_obj):
        """
        Publish dataset to be available for client figure to plot

        Parameters
        ----------
        name: :obj:`str`
            Name of the dataset
        data_source_obj: :class:`~pyhdx.panel.data_sources.DataSource`
            Data source object
        """

        try:  # update existing source
            src = self.sources[name]
            #todo next callback??
            src.source.data.update(**data_source_obj.source.data)  #todo refactor source to cds?
        except KeyError:
            self.sources[name] = data_source_obj

        self.param.trigger('sources')


class PyHDXController(MainController):
    """
    Main controller for PyHDX web application.

    """
    fit_results = param.Dict({}, doc='Dictionary of fit results', precedence=-1)
    peptides = param.ClassSelector(PeptideMasterTable, doc='Master list of all peptides', precedence=-1)
    series = param.ClassSelector(KineticsSeries,
                                 doc='KineticsSeries object with current selected and corrected peptides', precedence=-1)

    def __init__(self, *args, **kwargs):
        super(PyHDXController, self).__init__(*args, **kwargs)


class ComparisonController(MainController):
    """
    Main controller for binary comparison web application.
    """

    datasets = param.Dict(default={}, doc='Dictionary for all datasets')
    comparisons = param.Dict(default={}, doc='Dictionary for all comparisons (should be in sources)')
