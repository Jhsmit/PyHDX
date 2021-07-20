from pyhdx.fileIO import csv_to_dataframe
from pathlib import Path
from pyhdx.web.views import hvRectangleAppView
from pyhdx.web.transforms import PeptideLayoutTransform
from pyhdx.web.filters import UniqueValuesFilter, MultiIndexSelectFilter
from pyhdx.web.sources import DataFrameSource
from pyhdx.web.opts import CmapOpts
from bokeh.models.tools import HoverTool
import panel as pn

pn.extension(sizing_mode='stretch_both')


current_dir = Path(__file__).parent
data_dir = current_dir / 'test_data'

df_peptides = csv_to_dataframe(data_dir / 'peptides.txt')

tables = {'peptides': df_peptides}
source = DataFrameSource(tables=tables, name='dataframe')

hover = HoverTool(tooltips=[("index", "@index"), ('rfu', '@value (@uptake_corrected D)'),
                            ('sequence', '@sequence')])
additional_opts = {'color': 'value', 'colorbar': True, 'responsive': True, 'clim': (0, 1), 'framewise': True,
                   'xlabel': "Residue Number", 'ylabel': '', 'yticks': 0, 'tools': [hover]}
cmap_opts = CmapOpts(opts=additional_opts, name='cmap')

opts_list = [cmap_opts]

peptides_transform = PeptideLayoutTransform(
    value='rfu', name='trs_peptides',
    passthrough=['uptake', 'uptake_corrected', 'sequence', 'uptake_corrected', 'ex_residues'])

multiindex_select_filter = MultiIndexSelectFilter(field='state', name='select_index', table='peptides',
                                                  source=source)

slider_exposure_filter = UniqueValuesFilter(field='exposure', name='exposure_slider',
                                            table='peptides', filters=[multiindex_select_filter], source=source)
filters = [multiindex_select_filter, slider_exposure_filter]

coverage = hvRectangleAppView(source=source, name='coverage', table='peptides', opts=cmap_opts.opts,
                              streaming=True,
                              transforms=[peptides_transform],
                              filters=[multiindex_select_filter, slider_exposure_filter])

coverage.update()

def update(*events):
    coverage.update()

for f in filters:
    f.param.watch(update, 'value')

buttons = pn.Column(*[f.panel for f in filters])
row = pn.Row(*[buttons, coverage.panel])
pn.serve(row)