from pyhdx.panel.main_controllers import PyHDXController
from pyhdx.panel.controllers import *
from pyhdx.panel.base import BokehFigurePanel, STATIC_DIR
from pyhdx.panel.views import *
from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.log import get_default_handler
from pyhdx import VERSION_STRING_SHORT
from pyhdx.fileIO import csv_to_dataframe
from pyhdx.panel.sources import DataFrameSource
from pyhdx.panel.transforms import RescaleTransform, ApplyCmapTransform

import sys

import panel as pn
from panel import pane
from lumen.views import PerspectiveView, hvPlotView

from pathlib import Path
import matplotlib as mpl


"""
Example to test Lumen components to construct the web application

"""

current_dir = Path(__file__).parent
data_dir = current_dir.parent.parent / 'tests' / 'test_data'


control_panels = [
    CSVFileInputControl
]

df = csv_to_dataframe(data_dir / 'ecSecB_torch_fit.txt')

source = DataFrameSource(df=df, name='torch_fit')

print(source.get_schema())

rescale_transform = RescaleTransform(field='deltaG', scale_factor=1e-3)

cmap = mpl.cm.get_cmap('viridis')
norm = mpl.colors.Normalize(vmin=0, vmax=20)
cmap_transform = ApplyCmapTransform(cmap=cmap, norm=norm, field='deltaG')


# plot = hvPlotAppView(source=source, x='r_number', y='deltaG', kind='scatter', name='hvplot', c='color',
#                   table='torch_fit', transforms=[rescale_transform, cmap_transform])
#
# plot_panel = plot.get_panel()
# pn.serve(plot_panel)
