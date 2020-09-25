from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.controller import *
from pyhdx.panel.fig_panels import *
from pyhdx.panel.log import get_default_handler
import sys
from pyhdx import VERSION_STRING_SHORT
import panel as pn

DEBUG = True

control_panels = [
    PeptideFileInputControl,
    CoverageControl,
    InitialGuessControl,
    FitControl,
    FitResultControl,
    ClassificationControl,
    FileExportControl,
    ProteinViewControl,
    OptionsControl
]

if DEBUG:
    control_panels.append(DeveloperControl)

figure_panels = [
    CoverageFigure,
    RateFigure,
    PFactFigure,
    FitResultFigure,
    ProteinFigure,
    LoggingFigure
]

elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
cluster = '127.0.0.1:52123'
ctrl = PyHDXController(control_panels, figure_panels, cluster=cluster)
ctrl.logger.addHandler(get_default_handler(sys.stdout))
tmpl = elvis.compose(ctrl.control_panels.values(),
                     elvis.column(
                         elvis.stack(
                             elvis.view(ctrl.figure_panels['CoverageFigure']),
                             elvis.view(ctrl.figure_panels['ProteinFigure'])
                         ),
                         elvis.stack(
                             elvis.view(ctrl.figure_panels['RateFigure']),
                             elvis.view(ctrl.figure_panels['PFactFigure']),
                             elvis.view(ctrl.figure_panels['FitResultFigure']),
                             elvis.view(ctrl.figure_panels['LoggingFigure']),
                         )
                     )
                    )

ctrl.control_panels['OptionsControl']._update_link()

tmpl.servable(title='main')


