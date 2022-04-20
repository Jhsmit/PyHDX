"""
Reload SecB and fitted data and launch  GUI
Run local_cluster.py in anothor thread

"""

import sys
from pathlib import Path

import pandas as pd
import panel as pn
import yaml
import numpy as np

from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import csv_to_dataframe, load_fitresult
from pyhdx.fileIO import csv_to_protein
from pyhdx.web.apps import main_app
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.web.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.web.utils import load_state
from pyhdx.web.constructor import AppConstructor
from dev_controllers import *

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

cwd = Path(__file__).parent.resolve()
yaml_dict = yaml.safe_load((cwd / "mwe_app.yaml").read_text(encoding="utf-8"))

ctr = AppConstructor()
ctrl = ctr.parse(yaml_dict)

fmt = {"accent_base_color": "#1d417a"}

# tmpl = pn.template.FastGridTemplate(title=f'MWE app', **fmt)
# controllers = ctrl.control_panels.values()
# controls = pn.Accordion(*[controller.panel for controller in controllers], toggle=True)
# tmpl.sidebar.append(controls)

views_names = ["xy_scatter", "xy_line"]

elvis = GoldenElvis(
    ctrl, ExtendedGoldenTemplate, ExtendedGoldenDefaultTheme, title="test"
)


views = {v: ctrl.views[v] for v in views_names}
# [v.update() for v in views.values()]

tmpl = elvis.compose(elvis.row(elvis.view("xy_scatter"), elvis.view("xy_line")))

# cov_tab = pn.Tabs(
#     ('Coverage', views['coverage'].panel),
#     ('Protein', views['protein'].panel)
# )
#
# scatter_tab = pn.Tabs(
#     ('RFU', views['rfu_scatter'].panel),
#     ('Rates', views['rates'].panel),
#     ('ΔG', views['gibbs'].panel),
#     ('ΔΔG', views['ddG'].panel),
# )
#
# log_tab = pn.Tabs(
#     ('Info log', views['logging_info'].panel),
#     ('Debug log', views['logging_debug'].panel)
# )

# tmpl.main[0:3, 0:6] = views['xy_scatter'].panel
# tmpl.main[0:3, 6:12] = views['xy_line'].panel


def reload_tables():

    df = pd.DataFrame(
        {
            "x": np.random.normal(loc=3, scale=2, size=100),
            "y": np.random.normal(loc=2, scale=0.3, size=100),
        }
    )

    src = ctrl.sources["main"]
    src.tables["test_data"] = df
    src.param.trigger("updated")


# update tab sizes:
# https://stackoverflow.com/questions/39534178/golden-layout-how-to-increase-header-tabs-height/40390999
def reload_dashboard():
    pass


def init_dashboard():
    pass


# pn.state.onload(reload_dashboard)
pn.state.onload(reload_tables)
# pn.state.onload(init_dashboard)

if __name__ == "__main__":
    pn.serve(tmpl, show=True, static_dirs={"pyhdx": STATIC_DIR})

elif __name__.startswith("bokeh_app"):
    tmpl.servable()

# ctrl.template.servable()


# panel serve --show --autoreload --static-dirs pyhdx=C:\Users\jhsmi\pp\PyHDX\pyhdx\web\static
