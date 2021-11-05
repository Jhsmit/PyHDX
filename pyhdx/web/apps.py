from pyhdx.web.controllers import *
import sys
from pathlib import Path

import panel as pn
import yaml

from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import dataframe_to_file, csv_to_dataframe, load_fitresult
from pyhdx.local_cluster import default_client
from pyhdx.models import HDXMeasurementSet
from pyhdx.plot import default_cmap_norm
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.controllers import *
from pyhdx.web.filters import CrossSectionFilter, RescaleFilter, ApplyCmapOptFilter, TableSourceFilter
from pyhdx.web.log import get_default_handler
from pyhdx.web.log import logger
from pyhdx.web.main_controllers import PyHDXController
from pyhdx.web.opts import CmapOpts
from pyhdx.web.sources import PyHDXSource
from pyhdx.web.views import *
from pyhdx.web.views import hvRectanglesAppView

DEBUG = False

regen_data = True

cwd = Path(__file__).parent
root = cwd.parent.parent
data_input_dir = root / 'tests' / 'test_data' / 'input'
data_output_dir = root / 'tests' / 'test_data' / 'output'

print(data_output_dir)
yaml_stream = Path(data_input_dir / 'data_states.yaml').read_text()
data_dict = yaml.safe_load(yaml_stream)

test_dir = cwd / 'new_test_data'
test_dir.mkdir(exist_ok=True)

SCATTER_KWARGS = {'size': 8}

pn.extension(loading_spinner='dots', loading_color='#00aa41')

if regen_data:
    src = PyHDXSource(name='PyHDX')

    state = 'SecB_tetramer'
    hdxm_tetramer = yaml_to_hdxm(data_dict[state], data_dir=data_input_dir, name=state)
    src.add(hdxm_tetramer, state)

    state = 'SecB_dimer'
    hdxm_dimer = yaml_to_hdxm(data_dict[state], data_dir=data_input_dir, name=state)
    src.add(hdxm_dimer, state)
    dataframe_to_file(test_dir / 'peptides.csv', src.tables['peptides'])
    dataframe_to_file(test_dir / 'rfu_residues.csv', src.tables['rfu_residues'])

    hdxm_set = HDXMeasurementSet([hdxm_tetramer, hdxm_dimer])
    results = [fit_rates_half_time_interpolate(hdxm) for hdxm in hdxm_set]
    guess_result = RatesFitResult(results)
    src.add(guess_result, 'Guess_1')
    dataframe_to_file(test_dir / 'rates.csv', src.tables['rates'])

    fit_result = load_fitresult(data_output_dir / 'ecsecb_tetramer_dimer')
    src.add(fit_result, 'Fit_1')
    dataframe_to_file(test_dir / 'dG_fit.csv', src.tables['dG_fits'])


global_opts = {'show_grid': True}


@logger('pyhdx')
def main_app(client='default'):
    client = default_client() if client == 'default' else client
    logger = main_app.logger

    sources, transforms, filters, views, opts = {}, {}, {}, {}, {}


    src = PyHDXSource(name='PyHDX')
    # preload data

    src.tables['peptides'] = csv_to_dataframe(test_dir / 'peptides.csv')
    src.tables['rfu_residues'] = csv_to_dataframe(test_dir / 'rfu_residues.csv')
    src.tables['dG_fits'] = csv_to_dataframe(test_dir / 'dG_fit.csv')
    src.tables['rates'] = csv_to_dataframe(test_dir / 'rates.csv')

    sources[src.name] = src

    pdb_string = (test_dir / '1qyn.pdb').read_text()


    # ---------------------------------------------------------------------- #
    #                                Coverage Figure
    # ---------------------------------------------------------------------- #

    filters = {}

    f = TableSourceFilter(source=src, table='peptides', name='peptide_src')
    filters[f.name] = f

    f = CrossSectionFilter(n_levels=2, source=filters['peptide_src'],
                           name='coverage_select')
    filters[f.name] = f

    hover = HoverTool(tooltips=[("id", "@index"), ('rfu', '@rfu (@uptake_corrected D)'),
                                ('sequence', '@sequence')])
    static_opts = {'color': 'rfu', 'colorbar': True, 'responsive': True, 'framewise': True,  # what are framewise and responsive? stream / resize?
                   'xlabel': "Residue Number", 'ylabel': '', 'yticks': 0, 'tools': [hover], **global_opts}
    o = CmapOpts(name='rfu', field='rfu')
    opts[o.name] = o

    v = hvRectanglesAppView(
        source=filters['coverage_select'], name='coverage',
        opts={**opts['rfu'].opts, **static_opts},
        left='start', right='end', passthrough=['rfu', 'uptake_corrected', 'sequence'])
    views[v.name] = v


    # ---------------------------------------------------------------------- #
    #                                rates figure
    # ---------------------------------------------------------------------- #
    f = TableSourceFilter(source=src, table='rates', name='rates_src')
    filters[f.name] = f

    f = CrossSectionFilter(n_levels=2, source=filters['rates_src'],
                           name='rates_select')
    filters[f.name] = f

    static_opts = dict(logy=True, ylabel='Rate (s⁻¹)',
                 xlabel='Residue Number', invert_yaxis=True,
                       responsive=True,
                       **SCATTER_KWARGS,
                       **global_opts)

    v = hvScatterAppView(
        source=filters['rates_select'], name='rates', y='rate', x='r_number',
        opts=static_opts,
        ) #issue 154: deltaG units (Shady)
    views[v.name] = v

    # ---------------------------------------------------------------------- #
    #                                deltaG figure
    # ---------------------------------------------------------------------- #

    f = TableSourceFilter(source=src, table='dG_fits', name='dG_src')
    filters[f.name] = f

    f = CrossSectionFilter(n_levels=2, source=filters['dG_src'], name='dG_fit_select')
    filters[f.name] = f

    f = RescaleFilter(column='deltaG', scale_factor=1e-3, name='dG_rescale', source=filters['dG_fit_select'])
    filters[f.name] = f

    cmap, norm = default_cmap_norm('dG')
    norm.vmin *= 1e-3
    norm.vmax *= 1e-3
    static_opts = dict(color='deltaG', ylabel='ΔG (kJ/mol)',
                       xlabel='Residue Number',
                       invert_yaxis=True,
                       colorbar=True,
                       responsive=True,
                       **global_opts,
                       **SCATTER_KWARGS)
    o = CmapOpts(name='dG', cmap=cmap, norm=norm, field='deltaG', sclf=1e3)

    opts[o.name] = o

    v = hvScatterAppView(
        source=filters['dG_rescale'], name='gibbs', y='deltaG', x='r_number',
        opts={**opts['dG'].opts, **static_opts},
        )
    views[v.name] = v

    # ---------------------------------------------------------------------- #
    #                                Protein view
    # ---------------------------------------------------------------------- #

    f = TableSourceFilter(source=src, name='protein_table', table_options=['rfu_residues', 'dG_fits'])  # todo metaclass autoregister filters?
    filters[f.name] = f

    f = CrossSectionFilter(source=filters['protein_table'], name='protein_select', n_levels=0)
    filters[f.name] = f

    cmap_opts = {k: v for k, v in opts.items() if isinstance(v, CmapOpts)}
    f = ApplyCmapOptFilter(source=filters['protein_select'],
                           opts_dict = cmap_opts,
                           name='protein_cmap')  #todo allow concat filter which combined multiple filters
    filters[f.name] = f

    v = NGLView(
        source=filters['protein_cmap'],  # update to what it should actually be
        name='protein',
        object=pdb_string,
        dependencies=list(cmap_opts.values()),
    )
    views[v.name] = v

    # ---------------------------------------------------------------------- #
    #                                VIEWS
    # ---------------------------------------------------------------------- #

    control_panels = [
        DevTestControl,
        PeptideFileInputControl,
        InitialGuessControl,
        FitControl,
        ColorTransformControl,
        ProteinControl,
        FileExportControl,
        FigureExportControl,

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

    ctrl.src.param.trigger('updated')


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

    tmpl = pn.template.FastGridTemplate()
    controllers = ctrl.control_panels.values()
    controls = pn.Accordion(*[controller.panel for controller in controllers], toggle=True)
    tmpl.sidebar.append(controls)

    views_names = [
        'protein',
        'coverage',
        'rates',
        'gibbs'
        ]

    added_views = [ctrl.views[v] for v in views_names]

    for i, view in enumerate(added_views):
        view.update()
        item = pn.Row(view.panel, sizing_mode='stretch_both')
        tmpl.main[3*i:3*(i+1), :] = item

    pn.serve(tmpl, static_dirs={'pyhdx': STATIC_DIR})

#    pn.serve(ctrl.template, static_dirs={'pyhdx': STATIC_DIR})