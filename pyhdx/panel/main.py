from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.controller import Controller
from pyhdx import VERSION_STRING_SHORT
import panel as pn


control_panels = [
    'FileInputControl',
    'CoverageControl',
    'FittingControl',
    'TFFitControl',
    'FittingQuality',
    'ClassificationControl',
    'FileExportPanel',
    'OptionsPanel'
]

figure_panels = [
    'CoverageFigure',
    'RateFigure',
    'PFactFigure',
    'FitResultFigure'
]

elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
cluster = '127.0.0.1:52123'
ctrl = Controller(control_panels, figure_panels, cluster=cluster)

tmpl = elvis.compose(ctrl.control_panels.values(),
                              elvis.column(
                                  elvis.view(ctrl.figure_panels['CoverageFigure']),
                                  elvis.stack(
                                      elvis.view(ctrl.figure_panels['RateFigure']),
                                      elvis.view(ctrl.figure_panels['PFactFigure']),
                                      elvis.view(ctrl.figure_panels['FitResultFigure'])
                                  )
                              ))


if __name__ == '__main__':
    pn.serve(tmpl)


