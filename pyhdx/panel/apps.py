from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.controllers import *
from pyhdx.panel.main_controllers import ComparisonController, PyHDXController
from pyhdx.panel.views import *
from pyhdx.panel.log import get_default_handler
import sys
from pyhdx import VERSION_STRING_SHORT
from pyhdx.panel.base import BokehFigurePanel, STATIC_DIR
from pyhdx.fileIO import csv_to_dataframe
from pyhdx.panel.sources import DataFrameSource
from pyhdx.panel.transforms import RescaleTransform, ApplyCmapTransform, PeptideLayoutTransform, ResetIndexTransform
from pyhdx.panel.opts import CmapOpts
from pyhdx.panel.filters import UniqueValuesFilter, MultiIndexSelectFilter
from pyhdx.panel.log import StreamToLogger
import logging
import panel as pn
from pyhdx.panel.log import logger
from pyhdx.panel.config import ConfigurationSettings
from panel import pane
from lumen.views import PerspectiveView, hvPlotView
from lumen.filters import WidgetFilter, ParamFilter

from pathlib import Path
import pandas as pd
import matplotlib as mpl
import datetime

DEBUG = True

current_dir = Path(__file__).parent
data_dir = current_dir.parent.parent / 'tests' / 'test_data'
global_opts = {'show_grid': True}
cfg = ConfigurationSettings()

@logger('pyhdx')
def main_app():
    logger = main_app.logger

    # ---------------------------------------------------------------------- #
    #                                SOURCES
    # ---------------------------------------------------------------------- #

    col_index = pd.MultiIndex.from_tuples([], names=('state', 'quantity'))
    df_peptides = pd.DataFrame(columns=col_index)

    col_index = pd.MultiIndex.from_tuples([], names=('fit ID', 'state', 'quantity'))
    row_index = pd.RangeIndex(0, 1, name='r_number')
    df_rates = pd.DataFrame(columns=col_index, index=row_index)
    # todo make sure that proper-shaped df is used to initiate stream (and generalize for rectangles plot)

    col_index = pd.MultiIndex.from_tuples([], names=('fit ID', 'state', 'quantity'))
    row_index = pd.RangeIndex(0, 1, name='r_number')
    df_global_fit = pd.DataFrame(columns=col_index, index=row_index)

    # Availble tables are predefined at launch, but are empty
    # this way GUI methods can add to them as multiindex subset
    # more tables can be added later by the gui
    tables = {'peptides': df_peptides, 'rates': df_rates, 'global_fit': df_global_fit}
    source = DataFrameSource(tables=tables, name='dataframe')

    #df = csv_to_dataframe(data_dir / 'ecSecB_apo_peptides.csv')
    #source.add_df(df, 'peptides', 'ecSecB_apo')

    src_list = [source]

    # ---------------------------------------------------------------------- #
    #                                TRANSFORMS
    # ---------------------------------------------------------------------- #

    # rescale_transform = RescaleTransform(field='deltaG', scale_factor=1e-3)
    #
    # cmap = mpl.cm.get_cmap('viridis')
    # norm = mpl.colors.Normalize(vmin=0, vmax=20)
    # cmap_transform = ApplyCmapTransform(cmap=cmap, norm=norm, field='deltaG')

    peptides_transform = PeptideLayoutTransform(value='scores')

    reset_index_transform = ResetIndexTransform()

    trs_list = [peptides_transform, reset_index_transform]

    # ---------------------------------------------------------------------- #
    #                                FILTERS
    # ---------------------------------------------------------------------- #

    #unique_vals = list(np.sort(peptides_source.get_unique(table='peptides', field='exposure')))
    multiindex_select_filter = MultiIndexSelectFilter(field='state', name='select_index', table='peptides',
                                                      source=source)

     # unique_vals = list(np.sort(source.get_unique(table='peptides', field='exposure', state='ecSecB_apo')))
    slider_exposure_filter = UniqueValuesFilter(field='exposure', name='exposure_slider',
                                                table='peptides', filters=[multiindex_select_filter], source=source)

    filter_list = [multiindex_select_filter, slider_exposure_filter]

    # ---------------------------------------------------------------------- #
    #                                OPTS
    # ---------------------------------------------------------------------- #

    additional_opts = {'color': 'value', 'colorbar': True, 'responsive': True, 'clim': (0, 100), 'framewise': True,
                       'xlabel': "Residue Number", 'ylabel': '', 'yticks': 0, **global_opts}
    cmap_opts = CmapOpts(opts=additional_opts, name='cmap')

    opts_list = [cmap_opts]

    # ---------------------------------------------------------------------- #
    #                                VIEWS
    # ---------------------------------------------------------------------- #

    view_list = []

    rescale_transform = RescaleTransform(field='deltaG', scale_factor=1e-3)

    cmap = mpl.cm.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=0, vmax=20)
    cmap_transform = ApplyCmapTransform(cmap=cmap, norm=norm, field='deltaG', name='cmap_transform')

    trs_list.append(rescale_transform)
    trs_list.append(cmap_transform)

    multiindex_select_global_fit_1 = MultiIndexSelectFilter(field='fit ID', name='select_index_global_fit_lv1', table='global_fit',
                                                       source=source)
    multiindex_select_global_fit_2 = MultiIndexSelectFilter(field='state', name='select_index_global_fit_lv2', table='global_fit',
                                                       source=source, filters=[multiindex_select_global_fit_1])

    filter_list += [multiindex_select_global_fit_1, multiindex_select_global_fit_2]

    opts = {'xlabel': 'Residue Number', 'ylabel': 'ΔG (kJ mol⁻¹)', **global_opts}
    deltaG = hvPlotAppView(source=source, name='gibbs', x='r_number', y='deltaG', kind='scatter', c='color',
                           table='global_fit', transforms=[cmap_transform, rescale_transform], streaming=True,
                           filters = [multiindex_select_global_fit_1, multiindex_select_global_fit_2],
                           responsive=True, opts=opts) #issue 154: deltaG units


    view_list.append(deltaG)

    coverage = hvRectangleAppView(source=source, name='coverage', table='peptides', opts=cmap_opts.opts,
                                  streaming=True,
                                  transforms=[peptides_transform],
                                  filters=[multiindex_select_filter, slider_exposure_filter])
    view_list.append(coverage)



    multiindex_select_rates_1 = MultiIndexSelectFilter(field='fit ID', name='select_index_rates_lv1', table='rates',
                                                       source=source)
    multiindex_select_rates_2 = MultiIndexSelectFilter(field='state', name='select_index_rates_lv2', table='rates',
                                                       source=source, filters=[multiindex_select_rates_1])

    filter_list += [multiindex_select_rates_1, multiindex_select_rates_2]
    # perhaps consider derivedsource for the views

    opts = {'logy': True, 'xlabel': "Residue Number", 'ylabel': "Rate (min⁻¹)", **global_opts}
    rates = hvPlotAppView(source=source, name='rates', x='r_number', y='rate', kind='scatter', # c='color'
                           table='rates', streaming=True, responsive=True, opts=opts,
                          transforms=[reset_index_transform],
                          filters=[multiindex_select_rates_1, multiindex_select_rates_2]
                           )
    view_list.append(rates)

    log_view = LoggingView(logger=logger, level=logging.INFO, name='Info log')
    debug_log_view = LoggingView(logger=logger, level=logging.DEBUG, name='Debug log')
    view_list.append(log_view)
    view_list.append(debug_log_view)

    sources = {src.name: src for src in src_list}
    transforms = {trs.name: trs for trs in trs_list}
    filters = {filt.name: filt for filt in filter_list}
    views = {view.name: view for view in view_list}
    opts = {opt.name: opt for opt in opts_list}

    control_panels = [
        PeptideFileInputControl,
        CoverageControl,
        InitialGuessControl,
        FitControl,
        ClassificationControl,
        GraphControl,
        # FitResultControl,
        FileExportControl,
        # ProteinViewControl,
        # OptionsControl
    ]

    if DEBUG:
        control_panels.append(DeveloperControl)

    ctrl = PyHDXController(control_panels,
                           sources=sources,
                           transforms=transforms,
                           filters=filters,
                           opts=opts,
                           views=views,
                           cluster=cfg.cluster,
                           logger=logger)

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme,
                        title=VERSION_STRING_SHORT)

    ctrl.logger.addHandler(get_default_handler(sys.stdout))

    elvis.compose(ctrl, elvis.column(
        elvis.view(ctrl.views['coverage']),
        elvis.row(
            elvis.stack(
                elvis.view(ctrl.views['rates']),
                elvis.view(ctrl.views['gibbs']),
            ),
            elvis.stack(
                elvis.view(ctrl.views['Info log']),
                elvis.view(ctrl.views['Debug log'])
            )
        )
        )
    )

    return ctrl

#
# def single_app():
#     control_panels = [
#         SingleMappingFileInputControl,
#         ClassificationControl,
#         ProteinViewControl,
#         DifferenceFileExportControl,
#         OptionsControl,
#         DeveloperControl
#     ]
#
#     if DEBUG:
#         control_panels.append(DeveloperControl)
#
#     figure_panels = [
#         BinaryComparisonFigure,
#         ProteinFigure,
#         LoggingFigure
#     ]
#
#     elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
#     ctrl = ComparisonController(control_panels, figure_panels, cluster=cluster)
#     ctrl.logger.addHandler(get_default_handler(sys.stdout))
#     elvis.compose(ctrl,
#                   elvis.column(
#                       elvis.stack(
#                           elvis.view(ctrl.figure_panels['ProteinFigure'])
#                       ),
#                       elvis.row(
#                           elvis.stack(
#                              elvis.view(ctrl.figure_panels['BinaryComparisonFigure']),
#                           ),
#                           elvis.view(ctrl.figure_panels['LoggingFigure']),
#                       )
#                   ))
#
#     ctrl.control_panels['ClassificationControl'].log_space = False
#
#     return ctrl
#
#
# def diff_app():
#     control_panels = [
#         MappingFileInputControl,
#         DifferenceControl,
#         ClassificationControl,
#         ProteinViewControl,
#         DifferenceFileExportControl,
#         OptionsControl,
#         DeveloperControl
#     ]
#
#     if DEBUG:
#         control_panels.append(DeveloperControl)
#
#     figure_panels = [
#         BinaryComparisonFigure,
#         SingleValueFigure,
#         ProteinFigure,
#         LoggingFigure
#     ]
#
#     elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
#     ctrl = ComparisonController(control_panels, figure_panels, cluster=cluster)
#     ctrl.logger.addHandler(get_default_handler(sys.stdout))
#     elvis.compose(ctrl,
#                   elvis.column(
#                       elvis.stack(
#                           elvis.view(ctrl.figure_panels['ProteinFigure'])
#                       ),
#                       elvis.row(
#                           elvis.stack(
#                              elvis.view(ctrl.figure_panels['BinaryComparisonFigure']),
#                              elvis.view(ctrl.figure_panels['SingleValueFigure'])
#                           ),
#                           elvis.view(ctrl.figure_panels['LoggingFigure']),
#                       )
#                   ))
#
#     ctrl.control_panels['ClassificationControl'].log_space = False
#     return ctrl
#
#
# def folding_app():
#     control_panels = [
#         PeptideFoldingFileInputControl,
#         CoverageControl,
#         FoldingFitting,
#         FitResultControl,
#         ClassificationControl,
#         ProteinViewControl,
#         FileExportControl,
#         OptionsControl
#     ]
#
#     if DEBUG:
#         control_panels.append(DeveloperControl)
#
#     figure_panels = [
#         CoverageFigure,
#         RateFigure,
#         ScoresFigure,
#         FitResultFigure,
#         ProteinFigure,
#         LoggingFigure
#     ]
#
#     elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
#     ctrl = PyHDXController(control_panels, figure_panels, cluster=cluster)
#     ctrl.logger.addHandler(get_default_handler(sys.stdout))
#     elvis.compose(ctrl,
#                   elvis.column(
#                       elvis.stack(
#                           elvis.view(ctrl.figure_panels['CoverageFigure']),
#                           elvis.view(ctrl.figure_panels['ProteinFigure'])
#                       ),
#                       elvis.row(
#                           elvis.stack(
#                               elvis.view(ctrl.figure_panels['RateFigure']),
#                               elvis.view(ctrl.figure_panels['ScoresFigure']),
#                               elvis.view(ctrl.figure_panels['FitResultFigure'])
#                           ),
#                           elvis.view(ctrl.figure_panels['LoggingFigure']),
#                      )
#                   )
#                   )
#
#     ctrl.control_panels['ClassificationControl'].log_space = False
#     return ctrl
#
#
# def full_deuteration_app():
#     control_panels = [
#         FDPeptideFileInputControl,
#         FDCoverageControl,
#         OptionsControl
#     ]
#
#     if DEBUG:
#         control_panels.append(DeveloperControl)
#
#     figure_panels = [
#         CoverageFigure,
#         LoggingFigure
#     ]
#
#     elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
#     ctrl = PyHDXController(control_panels, figure_panels, cluster=cluster)
#     ctrl.logger.addHandler(get_default_handler(sys.stdout))
#     elvis.compose(ctrl,
#                   elvis.column(
#                           elvis.view(ctrl.figure_panels['CoverageFigure']),
#                           elvis.view(ctrl.figure_panels['LoggingFigure']),
#                       ))
#
#     return ctrl
#
#
# def color_matrix_app():
#     control_panels = [
#         MatrixMappingFileInputControl,
#         ColoringControl,
#         ProteinViewControl,
#     ]
#
#     if DEBUG:
#         control_panels.append(DeveloperControl)
#
#     figure_panels = [
#         ImageFigure,
#         ProteinFigure,
#         LoggingFigure
#     ]
#
#     elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
#     ctrl = ComparisonController(control_panels, figure_panels, cluster=cluster)
#     ctrl.logger.addHandler(get_default_handler(sys.stdout))
#     elvis.compose(ctrl,
#                   elvis.column(
#                       elvis.stack(
#                           elvis.view(ctrl.figure_panels['ProteinFigure'])
#                       ),
#                       elvis.row(
#                           elvis.stack(
#                              elvis.view(ctrl.figure_panels['ImageFigure']),
#                           ),
#                           elvis.view(ctrl.figure_panels['LoggingFigure']),
#                       )
#                   ))
#
#     return ctrl


if __name__ == '__main__':
    ctrl = main_app()
    pn.serve(ctrl.template, static_dirs={'pyhdx': STATIC_DIR})