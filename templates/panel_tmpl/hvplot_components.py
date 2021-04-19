from pyhdx.panel.main_controllers import PyHDXController
from pyhdx.panel.controllers import CSVFileInputControl
from pyhdx.panel.base import BokehFigurePanel, STATIC_DIR
from pyhdx.panel.views import hvPlotAppView
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

current_dir = Path(__file__).parent
data_dir = current_dir.parent.parent / 'tests' / 'test_data'


"""
Example to test Lumen components to construct the web application

"""


control_panels = [
    CSVFileInputControl
]

df = csv_to_dataframe(data_dir / 'ecSecB_torch_fit.txt')


source = DataFrameSource(df=df, name='torch_fit')

rescale_transform = RescaleTransform(field='deltaG', scale_factor=1e-3)

cmap = mpl.cm.get_cmap('viridis')
norm = mpl.colors.Normalize(vmin=0, vmax=20)
cmap_transform = ApplyCmapTransform(cmap=cmap, norm=norm, field='deltaG')

plot = hvPlotAppView(source=source, x='r_number', y='deltaG', kind='scatter', name='hvplot', c='color',
                  table='torch_fit', transforms=[rescale_transform, cmap_transform])


table = source.get('torch_fit')
rescale = rescale_transform.apply(table)
trs_table = cmap_transform.apply(rescale)

hvplot = trs_table.hvplot(x='r_number', y='deltaG', kind='scatter')
import holoviews as hv
rectangles = hv.Rectangles([(0, 0, 1, 1), (2, 3, 4, 6), (0.5, 2, 1.5, 4), (2, 1, 3.5, 2.5)])
plot_panel = pn.pane.HoloViews(object=hvplot)

pn.serve(plot_panel)

