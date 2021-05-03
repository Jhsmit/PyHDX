from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.controllers import *
from pyhdx.panel.main_controllers import ComparisonController, PyHDXController
from pyhdx.panel.fig_panels import *
from pyhdx.panel.config import ConfigurationSettings
from pyhdx.panel.log import get_default_handler
import sys
from pyhdx import VERSION_STRING_SHORT
from pyhdx.panel.base import STATIC_DIR


DEBUG = True
cfg = ConfigurationSettings()


def main_app():
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
        DeltaGFigure,
        PFactFigure,
        ScoresFigure,
        FitResultFigure,
        ProteinFigure,
        LoggingFigure
    ]

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
    ctrl = PyHDXController(control_panels, figure_panels, cluster=cfg.cluster)
    ctrl.logger.addHandler(get_default_handler(sys.stdout))
    elvis.compose(ctrl,
                  elvis.column(
                      elvis.stack(
                          elvis.view(ctrl.figure_panels['CoverageFigure']),
                          elvis.view(ctrl.figure_panels['ProteinFigure'])
                      ),
                      elvis.stack(
                          elvis.view(ctrl.figure_panels['RateFigure']),
                          elvis.view(ctrl.figure_panels['DeltaGFigure']),
                          elvis.view(ctrl.figure_panels['PFactFigure']),
                          elvis.view(ctrl.figure_panels['FitResultFigure']),
                          elvis.view(ctrl.figure_panels['ScoresFigure']),
                          elvis.view(ctrl.figure_panels['LoggingFigure']),
                      )
                  ))

    ctrl.control_panels['OptionsControl']._update_link()
    return ctrl


def single_app():
    control_panels = [
        SingleMappingFileInputControl,
        ClassificationControl,
        ProteinViewControl,
        DifferenceFileExportControl,
        OptionsControl,
        DeveloperControl
    ]

    if DEBUG:
        control_panels.append(DeveloperControl)

    figure_panels = [
        BinaryComparisonFigure,
        ProteinFigure,
        LoggingFigure
    ]

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
    ctrl = ComparisonController(control_panels, figure_panels, cluster=cfg.cluster)
    ctrl.logger.addHandler(get_default_handler(sys.stdout))
    elvis.compose(ctrl,
                  elvis.column(
                      elvis.stack(
                          elvis.view(ctrl.figure_panels['ProteinFigure'])
                      ),
                      elvis.row(
                          elvis.stack(
                             elvis.view(ctrl.figure_panels['BinaryComparisonFigure']),
                          ),
                          elvis.view(ctrl.figure_panels['LoggingFigure']),
                      )
                  ))

    ctrl.control_panels['ClassificationControl'].log_space = False

    return ctrl


def diff_app():
    control_panels = [
        MappingFileInputControl,
        DifferenceControl,
        ClassificationControl,
        ProteinViewControl,
        DifferenceFileExportControl,
        OptionsControl,
        DeveloperControl
    ]

    if DEBUG:
        control_panels.append(DeveloperControl)

    figure_panels = [
        BinaryComparisonFigure,
        SingleValueFigure,
        ProteinFigure,
        LoggingFigure
    ]

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
    ctrl = ComparisonController(control_panels, figure_panels, cluster=cfg.cluster)
    ctrl.logger.addHandler(get_default_handler(sys.stdout))
    elvis.compose(ctrl,
                  elvis.column(
                      elvis.stack(
                          elvis.view(ctrl.figure_panels['ProteinFigure'])
                      ),
                      elvis.row(
                          elvis.stack(
                             elvis.view(ctrl.figure_panels['BinaryComparisonFigure']),
                             elvis.view(ctrl.figure_panels['SingleValueFigure'])
                          ),
                          elvis.view(ctrl.figure_panels['LoggingFigure']),
                      )
                  ))

    ctrl.control_panels['ClassificationControl'].log_space = False
    return ctrl


def folding_app():
    control_panels = [
        PeptideFoldingFileInputControl,
        CoverageControl,
        FoldingFitting,
        FitResultControl,
        ClassificationControl,
        ProteinViewControl,
        FileExportControl,
        OptionsControl
    ]

    if DEBUG:
        control_panels.append(DeveloperControl)

    figure_panels = [
        CoverageFigure,
        RateFigure,
        ScoresFigure,
        FitResultFigure,
        ProteinFigure,
        LoggingFigure
    ]

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
    ctrl = PyHDXController(control_panels, figure_panels, cluster=cfg.cluster)
    ctrl.logger.addHandler(get_default_handler(sys.stdout))
    elvis.compose(ctrl,
                  elvis.column(
                      elvis.stack(
                          elvis.view(ctrl.figure_panels['CoverageFigure']),
                          elvis.view(ctrl.figure_panels['ProteinFigure'])
                      ),
                      elvis.row(
                          elvis.stack(
                              elvis.view(ctrl.figure_panels['RateFigure']),
                              elvis.view(ctrl.figure_panels['ScoresFigure']),
                              elvis.view(ctrl.figure_panels['FitResultFigure'])
                          ),
                          elvis.view(ctrl.figure_panels['LoggingFigure']),
                     )
                  )
                  )

    ctrl.control_panels['ClassificationControl'].log_space = False
    return ctrl


def full_deuteration_app():
    control_panels = [
        FDPeptideFileInputControl,
        FDCoverageControl,
        OptionsControl
    ]

    if DEBUG:
        control_panels.append(DeveloperControl)

    figure_panels = [
        CoverageFigure,
        LoggingFigure
    ]

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
    ctrl = PyHDXController(control_panels, figure_panels, cluster=cfg.cluster)
    ctrl.logger.addHandler(get_default_handler(sys.stdout))
    elvis.compose(ctrl,
                  elvis.column(
                          elvis.view(ctrl.figure_panels['CoverageFigure']),
                          elvis.view(ctrl.figure_panels['LoggingFigure']),
                      ))

    return ctrl


def color_matrix_app():
    control_panels = [
        MatrixMappingFileInputControl,
        ColoringControl,
        ProteinViewControl,
    ]

    if DEBUG:
        control_panels.append(DeveloperControl)

    figure_panels = [
        ImageFigure,
        ProteinFigure,
        LoggingFigure
    ]

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
    ctrl = ComparisonController(control_panels, figure_panels, cluster=cfg.cluster)
    ctrl.logger.addHandler(get_default_handler(sys.stdout))
    elvis.compose(ctrl,
                  elvis.column(
                      elvis.stack(
                          elvis.view(ctrl.figure_panels['ProteinFigure'])
                      ),
                      elvis.row(
                          elvis.stack(
                             elvis.view(ctrl.figure_panels['ImageFigure']),
                          ),
                          elvis.view(ctrl.figure_panels['LoggingFigure']),
                      )
                  ))

    return ctrl


if __name__ == '__main__':
    ctrl = main_app()
    pn.serve(ctrl.template, static_dirs={'pyhdx': STATIC_DIR})
