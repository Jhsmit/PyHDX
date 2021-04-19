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

current_dir = Path(__file__).parent
data_dir = current_dir.parent.parent / 'tests' / 'test_data'


"""
Example to test Lumen components and application in Elvis template / web app

"""

def test_app():

    control_panels = [
        CSVFileInputControl,
        DeveloperControl
    ]

    df = csv_to_dataframe(data_dir / 'ecSecB_torch_fit.txt')

    source = DataFrameSource(df=df, name='torch_fit')

    rescale_transform = RescaleTransform(field='deltaG', scale_factor=1e-3)

    cmap = mpl.cm.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=0, vmax=20)
    cmap_transform = ApplyCmapTransform(cmap=cmap, norm=norm, field='deltaG')

    plot = hvPlotAppView(source=source, x='r_number', y='deltaG', kind='scatter', name='hvplot', c='color',
                         table='torch_fit', transforms=[rescale_transform, cmap_transform])

    figure_panels = [
        plot,
    ]

    cluster = 'none'

    elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
    ctrl = PyHDXController(control_panels, figure_panels, cluster=cluster)
    ctrl.logger.addHandler(get_default_handler(sys.stdout))

    elvis.compose(ctrl, elvis.column(
                  elvis.view(ctrl.figure_panels['hvplot']),
                  )

    )

    return ctrl

if __name__ == '__main__':

    ctrl = test_app()

    # root = Path(__file__).parent
    # test_file = root.parent / 'tests' / 'test_data' / 'ecSecB_torch_fit.txt'
    # print(test_file.read_text())
    #
    # binary = test_file.read_bytes()
    # input_control = ctrl.control_panels['CSVFileInputControl']
    # input_control.input_file = binary
    # input_control._action_load()

    #with open(root.parent / 'tests' / 'test_data' / 'ecSecB_torch_fit.txt'

    pn.serve(ctrl.template, static_dirs={'pyhdx': STATIC_DIR})