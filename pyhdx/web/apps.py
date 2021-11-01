from pyhdx.web.controllers import *
import sys
from pathlib import Path

import panel as pn

from pyhdx.fileIO import csv_to_dataframe
from pyhdx.local_cluster import default_client
from pyhdx.plot import default_cmap_norm
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.controllers import *
from pyhdx.web.filters import CrossSectionFilter, RescaleFilter
from pyhdx.web.log import get_default_handler
from pyhdx.web.log import logger
from pyhdx.web.main_controllers import PyHDXController
from pyhdx.web.opts import CmapOpts
from pyhdx.web.sources import PyHDXSource
from pyhdx.web.views import *
from pyhdx.web.views import hvRectangleAppView

DEBUG = True

current_dir = Path(__file__).parent
data_dir = current_dir.parent.parent / 'tests' / 'test_data'
global_opts = {'show_grid': True}

@logger('pyhdx')
def main_app(client='default'):
    client = default_client() if client == 'default' else client
    logger = main_app.logger

    sources, transforms, filters, views, opts = {}, {}, {}, {}, {}

    src = PyHDXSource(name='PyHDX')
    # preload data

    src.tables['peptides'] = csv_to_dataframe(test_dir / 'peptides.csv')
    src.tables['dG_fits'] = csv_to_dataframe(test_dir / 'dG_fit.csv')
    src.tables['rates'] = csv_to_dataframe(test_dir / 'rates.csv')

    sources[src.name] = src

    # ---------------------------------------------------------------------- #
    #                                Coverage Figure
    # ---------------------------------------------------------------------- #

    filters = {}

    #cs_filters = ['peptide_stack', 'peptide_pivot', 'peptide_reorder', 'peptide_sort_index']
    f = CrossSectionFilter(n_levels=2, source=src,
                           table='peptides', name='coverage_select')
    filters[f.name] = f

    hover = HoverTool(tooltips=[("id", "@index"), ('rfu', '@rfu (@uptake_corrected D)'),
                                ('sequence', '@sequence')])
    static_opts = {'color': 'rfu', 'colorbar': True, 'responsive': True, 'framewise': True,  # what are framewise and responsive? stream / resize?
                   'xlabel': "Residue Number", 'ylabel': '', 'yticks': 0, 'tools': [hover], **global_opts}
    o = CmapOpts(name='rfu')
    opts[o.name] = o

    v = hvRectangleAppView(
        source=src, name='coverage', table='peptides',
        opts={**opts['rfu'].opts, **static_opts},
        streaming=True,
        filters=[filters['coverage_select']],
        left='start', right='end', passthrough=['rfu', 'uptake_corrected', 'sequence'])
    views[v.name] = v


    # ---------------------------------------------------------------------- #
    #                                rates figure
    # ---------------------------------------------------------------------- #

    f = CrossSectionFilter(n_levels=2, source=src,
                           table='rates', name='rates_select')
    filters[f.name] = f

    static_opts = dict(logy=True, ylabel='Rate (s⁻¹)',
                 xlabel='Residue Number', invert_yaxis=True,
                       responsive=True,
                       **SCATTER_KWARGS)

    v = hvScatterAppView(
        source=src, name='rates', y='rate', x='r_number',
        table='rates', opts=static_opts,
        filters=[filters['rates_select']]) #issue 154: deltaG units (Shady)
    views[v.name] = v

    # ---------------------------------------------------------------------- #
    #                                deltaG figure
    # ---------------------------------------------------------------------- #

    f = CrossSectionFilter(n_levels=2, source=src, table='dG_fits', name='dG_fit_select')
    filters[f.name] = f

    f = RescaleFilter(column='deltaG', scale_factor=1e-3, name='dG_rescale')
    filters[f.name] = f

    cmap, norm = default_cmap_norm('dG')
    norm.vmin *= 1e-3
    norm.vmax *= 1e-3
    print(norm)
    static_opts = dict(color='deltaG', ylabel='ΔG (kJ/mol)',
                       xlabel='Residue Number',
                       invert_yaxis=True,
                       colorbar=True,
                       responsive=True,
                       **SCATTER_KWARGS)
    o = CmapOpts(name='dG', cmap=cmap, norm=norm)

    opts[o.name] = o

    v = hvScatterAppView(
        source=src, name='gibbs', y='deltaG', x='r_number',
        opts={**opts['dG'].opts, **static_opts},
        table='dG_fits', streaming=True, responsive=True,
        filters=[filters['dG_fit_select'], filters['dG_rescale']]) #issue 154: deltaG units (Shady)
    views[v.name] = v
    # view_list.append(deltaG)


    # ---------------------------------------------------------------------- #
    #                                VIEWS
    # ---------------------------------------------------------------------- #

    control_panels = [
        DevTestControl,
        PeptideFileInputControl,
        # CoverageControl,
        InitialGuessControl,
        FitControl,
        # GraphControl,
        ClassificationControl,
        # # ProteinControl,
        # # FitResultControl,
        # FileExportControl,
        # # ProteinViewControl,
        # # OptionsControlasdfasdf
    ]

    # if DEBUG:
    #     control_panels.append(DeveloperControl)

    ctrl = PyHDXController(control_panels,
                           src=src,
                           transforms=transforms,
                           filters=filters,
                           opts=opts,
                           views=views,
                           client=client,
                           logger=logger)

    ctrl.logger.addHandler(get_default_handler(sys.stdout))

    # elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme,
    #                     title=VERSION_STRING)
    #
    #
    # elvis.compose(ctrl,
    #               elvis.column(
    #     elvis.view(ctrl.views['coverage'])
    #               )
    #               )
    #
    return ctrl



if __name__ == '__main__':
    sys._excepthook = sys.excepthook

    import traceback as tb
    def my_exception_hook(exctype, value, traceback):
        # Print the error and traceback
        # https://stackoverflow.com/questions/43039048/pyqt5-fails-with-cryptic-message/43039363#43039363
        tb.print_tb(traceback, file=sys.stdout)
        print(exctype, value, traceback)

        tb.print_stack()
        print(traceback.format_exc())
        # or
        print(sys.exc_info()[2])
        # Call the normal Exception hook after
        sys._excepthook(exctype, value, traceback)
        sys.exit(1)


    # Set the exception hook to our wrapping function
    sys.excepthook = my_exception_hook

    ctrl = main_app()
    pn.serve(ctrl.template, static_dirs={'pyhdx': STATIC_DIR})