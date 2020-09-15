from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.controller import PyHDXController, ComparisonController
from pyhdx.panel.log import get_default_handler
import sys
from pyhdx import VERSION_STRING_SHORT
import panel as pn

DEBUG = False

control_panels = [
    'MappingFileInputControl',
    'DifferenceControl',
    'ClassificationControl',
    'FileExportControl',
    'ProteinViewControl',
    'OptionsControl'
]

if DEBUG:
    control_panels.append('DeveloperControl')

figure_panels = [
    #'LogLinearFigure',
    'PFactFigure',
    'ProteinFigure',
    'LoggingFigure'
]

elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
cluster = '127.0.0.1:52123'
ctrl = ComparisonController(control_panels, figure_panels, cluster=cluster)
ctrl.logger.addHandler(get_default_handler(sys.stdout))
tmpl = elvis.compose(ctrl.control_panels.values(),
                     elvis.column(
                         elvis.stack(
                             elvis.view(ctrl.figure_panels['ProteinFigure'])
                         ),
                         elvis.stack(
                             elvis.view(ctrl.figure_panels['PFactFigure']),
                             elvis.view(ctrl.figure_panels['LoggingFigure']),
                         )
                     )
                    )


if __name__ == '__main__':
    pn.serve(tmpl)


