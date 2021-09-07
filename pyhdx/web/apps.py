from pyhdx.web.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.web.theme import ExtendedGoldenDarkTheme
from pyhdx.web.controllers import *
from pyhdx.web.main_controllers import PyHDXController
from pyhdx.web.views import *
from pyhdx.web.log import get_default_handler
import sys
from pyhdx import VERSION_STRING
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.sources import DataFrameSource
from pyhdx.web.transforms import RescaleTransform, RemoveValueTransform, ApplyCmapTransform, PeptideLayoutTransform, ResetIndexTransform, \
    AccumulateRegularizersTransform
from pyhdx.web.opts import CmapOpts
from pyhdx.web.filters import UniqueValuesFilter, MultiIndexSelectFilter
import logging
import panel as pn
from pyhdx.web.log import logger
from pyhdx.config import ConfigurationSettings
from pyhdx.local_cluster import default_client

from pathlib import Path
import pandas as pd
import matplotlib as mpl

DEBUG = True

current_dir = Path(__file__).parent
data_dir = current_dir.parent.parent / 'tests' / 'test_data'
global_opts = {'show_grid': True}
cfg = ConfigurationSettings()

@logger('pyhdx')
def main_app(client='default'):
    client = default_client() if client == 'default' else client
    logger = main_app.logger

    # ---------------------------------------------------------------------- #
    #                                SOURCES
    # ---------------------------------------------------------------------- #

    tables = {}
    col_index = pd.MultiIndex.from_tuples([], names=('state_name', 'quantity'))
    df = pd.DataFrame(columns=col_index)
    tables['peptides'] = df

    col_index = pd.MultiIndex.from_tuples([], names=('fit_ID', 'state_name', 'quantity'))
    df = pd.DataFrame(columns=col_index)
    #df_peptides_mse.index.name = 'peptide index'
    tables['peptides_mse'] = df

    # Very annoying and hopefully temporary long-form dataframe
    col_index = pd.MultiIndex.from_tuples([], names=('fit_ID', 'state_name', 'quantity'))
    df = pd.DataFrame(columns=col_index)
    tables['d_calc'] = df

    col_index = pd.MultiIndex.from_tuples([], names=('state_name', 'exposure'))
    row_index = pd.RangeIndex(0, 1, name='r_number')
    df = pd.DataFrame(columns=col_index, index=row_index)
    tables['rfu'] = df

    col_index = pd.MultiIndex.from_tuples([], names=('fit_ID', 'state_name', 'quantity'))
    row_index = pd.RangeIndex(0, 1, name='r_number')
    df = pd.DataFrame(columns=col_index, index=row_index)
    tables['rates'] = df

    col_index = pd.MultiIndex.from_tuples([], names=('fit_ID', 'state_name', 'quantity'))
    row_index = pd.RangeIndex(0, 1, name='r_number')
    df = pd.DataFrame(columns=col_index, index=row_index)
    tables['global_fit'] = df

    col_index = pd.MultiIndex.from_tuples([], names=('fit_ID', 'state_name', 'loss_type'))
    row_index = pd.RangeIndex(0, 1, name='epochs')
    df = pd.DataFrame(columns=col_index, index=row_index)
    tables['losses'] = df

    col_index = pd.MultiIndex.from_tuples([], names=('color_ID', 'state_name', 'quantity'))
    row_index = pd.RangeIndex(0, 1, name='r_number')
    df = pd.DataFrame(columns=col_index, index=row_index)
    tables['colors'] = df

    # Availble tables are predefined at launch, but are empty
    # this way GUI methods can add to them as multiindex subset
    # tables = {'peptides': df_peptides,
    #           'peptides_mse': df_peptides_mse,
    #           'd_calc': df_d_calc,
    #           'rfu': df_rfu,
    #           'rates': df_rates,
    #           'global_fit': df_global_fit,
    #           'colors': df_colors}
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

    peptides_transform = PeptideLayoutTransform(
        value='rfu', name='trs_peptides',
        passthrough=['uptake', 'uptake_corrected', 'sequence', 'uptake_corrected', 'ex_residues', 'id']
    )
    peptides_mse_transform = PeptideLayoutTransform(value='total_mse', name='trs_peptides_mse')
    reset_index_transform = ResetIndexTransform(name='reset_index_trs')

    trs_list = [peptides_transform, peptides_mse_transform, reset_index_transform]


    # ---------------------------------------------------------------------- #
    #                                VIEWS
    # ---------------------------------------------------------------------- #

    view_list = []
    filters = {}
    # COVERAGE PEPTIDES VIEW
    f = MultiIndexSelectFilter(
        field='state_name', name='coverage_state_name', table='peptides', source=source)
    filters[f.name] = f

    f = UniqueValuesFilter(
        field='exposure', name='coverage_exposure', table='peptides', source=source,
        filters=[filters['coverage_state_name']])
    filters[f.name] = f

    hover = HoverTool(tooltips=[("index", "@id (@index)"), ('rfu', '@value (@uptake_corrected D)'),
                                ('sequence', '@sequence')])
    additional_opts = {'color': 'value', 'colorbar': True, 'responsive': True, 'clim': (0, 1), 'framewise': True,
                       'xlabel': "Residue Number", 'ylabel': '', 'yticks': 0, 'tools': [hover], **global_opts}
    cmap_opts = CmapOpts(opts=additional_opts, name='cmap')

    opts_list = [cmap_opts]

    coverage = hvRectangleAppView(
        source=source, name='coverage', table='peptides', opts=cmap_opts.opts, streaming=True,
        transforms=[peptides_transform],
        filters=[filters['coverage_state_name'], filters['coverage_exposure']])
    view_list.append(coverage)

    # # COVERAGE PEPTIDES MSE VIEW
    f = MultiIndexSelectFilter(
        field='fit_ID', name='coverage_mse_fit_id', table='peptides_mse', source=source)
    filters[f.name] = f

    f = MultiIndexSelectFilter(
        field='state_name', name='coverage_mse_state_name', table='peptides_mse', source=source,
        filters=[filters['coverage_mse_fit_id']])
    filters[f.name] = f

    hover = HoverTool(tooltips=[("index", "@index"), ('mse', '@value')])
    additional_opts_mse = {'color': 'value', 'colorbar': True, 'responsive': True, 'framewise': True,
                       'xlabel': "Residue Number", 'ylabel': '', 'yticks': 0, 'tools': [hover], **global_opts}
    cmap_mse_opts = CmapOpts(opts=additional_opts_mse, name='cmap_peptides_mse')
    coverage_mse = hvRectangleAppView(
        source=source, name='coverage_mse', table='peptides_mse', opts=cmap_mse_opts.opts, streaming=True,
        transforms=[peptides_mse_transform],
        filters=[filters['coverage_mse_fit_id'], filters['coverage_mse_state_name']])
    view_list.append(coverage_mse)

    # PEPTIDES VIEW vs time D_EXP
    logx = True
    opts = {'xlabel': 'Time (s)', 'ylabel': 'D-uptake', **global_opts}
    opts['logx'] = logx

    f = MultiIndexSelectFilter(
        field='state_name', name='peptide_d_exp_state_name', table='peptides', source=source)
    filters[f.name] = f

    f = UniqueValuesFilter(
        field='start_end', name='peptide_d_exp_select', show_index=True, table='peptides', source=source,
        filters=[filters['peptide_d_exp_state_name']])
    filters[f.name] = f

    removezero_trans = RemoveValueTransform(value=0., field='exposure')
    trs_list.append(removezero_trans)

    peptide_view = hvPlotAppView(
        source=source, name='peptide_view', table='peptides', streaming=True, kind='scatter',
        x='exposure', y='uptake_corrected', responsive=True, opts=opts,
        filters=[filters['peptide_d_exp_state_name'], filters['peptide_d_exp_select']],
        transforms=[removezero_trans])
    view_list.append(peptide_view)


    # PEPTIDES VIEW vs time D_CALC
    f = MultiIndexSelectFilter(
        field='fit_ID', name='peptide_d_calc_fit_id', table='d_calc', source=source)
    filters[f.name] = f

    f = MultiIndexSelectFilter(
        field='state_name', name='peptide_d_calc_state_name', table='d_calc', source=source,
        filters=[filters['peptide_d_calc_fit_id']])
    filters[f.name] = f

    f = UniqueValuesFilter(
        field='start_end', name='peptide_d_calc_select', show_index=True, table='d_calc', source=source,
        filters=[filters['peptide_d_calc_fit_id'], filters['peptide_d_calc_state_name']])
    filters[f.name] = f

    opts = {'xlabel': 'Time (s)', 'ylabel': 'D-uptake', 'color': 'r', **global_opts}
    opts['logx'] = logx
    d_calc_view = hvPlotAppView(
        source=source, name='d_calc_view', table='d_calc', streaming=True, kind='line',
        x='timepoints', y='d_calc', responsive=True, opts=opts,
        filters=[filters['peptide_d_calc_fit_id'], filters['peptide_d_calc_state_name'], filters['peptide_d_calc_select']])
    view_list.append(d_calc_view)

    # DELTA G VIEW
    rescale_transform = RescaleTransform(field='deltaG', scale_factor=1e-3)

    cmap = mpl.cm.get_cmap('viridis')
    cmap = mpl.colors.ListedColormap(["lawngreen"])
    norm = mpl.colors.Normalize(vmin=-1e6, vmax=1e6)
    cmap_transform = ApplyCmapTransform(cmap=cmap, norm=norm, field='deltaG', name='cmap_transform')

    trs_list.append(rescale_transform)
    trs_list.append(cmap_transform)

    f = MultiIndexSelectFilter(
        field='fit_ID', name='deltaG_fit_id', table='global_fit', source=source)
    filters[f.name] = f

    f = MultiIndexSelectFilter(
        field='state_name', name='deltaG_state_name', table='global_fit', source=source,
        filters=[filters['deltaG_fit_id']])
    filters[f.name] = f

    opts = {'xlabel': 'Residue Number', 'ylabel': 'ΔG (kJ mol⁻¹)', 'colorbar': False,  **global_opts}
    deltaG = hvPlotAppView(
        source=source, name='gibbs', x='r_number', y='deltaG', kind='scatter', c='color',
        table='global_fit', streaming=True, responsive=True, opts=opts,
        transforms=[cmap_transform, rescale_transform],
        filters=[filters['deltaG_fit_id'], filters['deltaG_state_name']],) #issue 154: deltaG units (Shady)

    view_list.append(deltaG)


    # RATES VIEW
    f = MultiIndexSelectFilter(
        field='fit_ID', name='rates_fit_id', table='rates', source=source)
    filters[f.name] = f

    f = MultiIndexSelectFilter(
        field='state_name', name='rates_state_name', table='rates', source=source,
        filters=[filters['rates_fit_id']])
    filters[f.name] = f

    # perhaps consider derivedsource for the views

    opts = {'logy': True, 'xlabel': "Residue Number", 'ylabel': "Rate (s⁻¹)", 'colorbar': False, **global_opts}
    rates = hvPlotAppView(
        source=source, name='rates', x='r_number', y='rate', kind='scatter', # c='color'
        table='rates', streaming=True, responsive=True, opts=opts,
        transforms=[reset_index_transform],
        filters=[filters['rates_fit_id'], filters['rates_state_name']]
                           )
    view_list.append(rates)


    # LOSSES VIEW

    f = MultiIndexSelectFilter(
        field='fit_ID', name='losses_fit_id', table='losses', source=source
    )
    filters[f.name] = f

    f = MultiIndexSelectFilter(
        field='state_name', name='losses_state_name', table='losses', source=source,
        filters=[filters['losses_fit_id']]
    )
    filters[f.name] = f

    reset_index_transform_loss = ResetIndexTransform(name='reset_index_losses')
    trs_list.append(reset_index_transform_loss)
    opts = {'xlabel': "Epochs", 'ylabel': "Loss", **global_opts}
    losses = hvPlotAppView(
        source=source, name='losses', x='index', y='mse_loss', kind='line', # c='color'
        table='losses', streaming=True, responsive=True, opts=opts, label='mse',
        transforms=[reset_index_transform_loss],
        filters=[filters['losses_fit_id'], filters['losses_state_name']]
    )
    view_list.append(losses)

    accumulate_reg_trnsform = AccumulateRegularizersTransform(name='accumulate_regularizers')
    trs_list.append(accumulate_reg_trnsform)

    opts = {'color': 'r', **opts}
    reg_losses = hvPlotAppView(
        source=source, name='reg_losses', x='index', y='reg_loss', kind='line',
        table='losses', streaming=True, responsive=True, opts=opts, label='reg',
        transforms=[reset_index_transform_loss, accumulate_reg_trnsform],
        filters=[filters['losses_fit_id'], filters['losses_state_name']]
    )
    view_list.append(reg_losses)


    # PROTEIN NGL VIEW
    f = MultiIndexSelectFilter(
        field='color_ID', name='ngl_color_id', table='colors', source=source
    )
    filters[f.name] = f

    f = MultiIndexSelectFilter(
        field='state_name', name='ngl_state_name', table='colors', source=source,
        filters=[filters['ngl_color_id']]
    )
    filters[f.name] = f

    protein_view = NGLView(
        source=source, name='protein', table='colors',
        filters=[filters['ngl_color_id'], filters['ngl_state_name']]
    )
    view_list.append(protein_view)

    # LOGGING VIEWS
    log_view = LoggingView(logger=logger, level=logging.INFO, name='Info log')
    debug_log_view = LoggingView(logger=logger, level=logging.DEBUG, name='Debug log')
    view_list.append(log_view)
    view_list.append(debug_log_view)

    sources = {src.name: src for src in src_list}
    transforms = {trs.name: trs for trs in trs_list}
    #filters = {filt.name: filt for filt in filter_list}
    views = {view.name: view for view in view_list}
    opts = {opt.name: opt for opt in opts_list}

    control_panels = [
        PeptideFileInputControl,
        CoverageControl,
        InitialGuessControl,
        FitControl,
        GraphControl,
        ClassificationControl,
        ProteinControl,
        # FitResultControl,
        FileExportControl,
        # ProteinViewControl,
        # OptionsControlasdfasdf
    ]

    if DEBUG:
        control_panels.append(DeveloperControl)

    ctrl = PyHDXController(control_panels,
                           sources=sources,
                           transforms=transforms,
                           filters=filters,
                           opts=opts,
                           views=views,
                           client=client,
                           logger=logger)

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme,
                        title=VERSION_STRING)

    ctrl.logger.addHandler(get_default_handler(sys.stdout))

    elvis.compose(ctrl, elvis.column(
        elvis.stack(
            elvis.view(ctrl.views['coverage']),
            elvis.view(ctrl.views['protein']),
            elvis.view(ctrl.views['coverage_mse'])
        ),
        elvis.row(
            elvis.stack(
                elvis.view(ctrl.views['rates']),
                elvis.view(ctrl.views['gibbs']),
                elvis.view([ctrl.views['peptide_view'], ctrl.views['d_calc_view']], title='peptide'),
                #elvis.view(ctrl.views['d_calc_view']),
            ),
            elvis.stack(
                elvis.view(ctrl.views['Info log']),
                elvis.view(ctrl.views['Debug log']),
                elvis.view([ctrl.views['losses'], ctrl.views['reg_losses']], title='losses')
                #elvis.view(ctrl.views['losses'])
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
    print('joe')

    sys._excepthook = sys.excepthook

    import traceback as tb
    def my_exception_hook(exctype, value, traceback):
        # Print the error and traceback

        tb.print_tb(traceback, file=sys.stdout)
        print(exctype, value, traceback)
        print('whooo')

        tb.print_stack()
        print(traceback.format_exc())

        print('whooo2')
        # or
        print(sys.exc_info()[2])
        # Call the normal Exception hook after
        sys._excepthook(exctype, value, traceback)
        sys.exit(1)


    # Set the exception hook to our wrapping function
    sys.excepthook = my_exception_hook

    ctrl = main_app()
    pn.serve(ctrl.template, static_dirs={'pyhdx': STATIC_DIR})