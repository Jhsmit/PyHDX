"""
Launch RFU app and preload data
Run local_cluster.py in anothor thread

"""
import sys
from pathlib import Path

import pandas as pd
import panel as pn
import yaml

from pyhdx.batch_processing import StateParser
from pyhdx.fileIO import csv_to_dataframe, load_fitresult
from pyhdx.fileIO import csv_to_protein
from pyhdx.web.apps import rfu_app
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.utils import load_state_rfu, fix_multiindex_dtypes
from pyhdx.config import cfg, reset_config

reset_config()

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


ctrl, tmpl = rfu_app()

cwd = Path(__file__).parent
root_dir = cwd.parent.parent
data_dir = root_dir / "tests" / "test_data" / "input"


# batch_fname = 'data_states.yaml' # standard secb apo / dimer dataset
# batch_fname = 'data_states_red.yaml'  # reduced number of pepties
batch_fname = "PpiX_states.yaml"  # secb apo / dimer but artificial delta C/N tail

state_spec = yaml.safe_load(Path(data_dir / batch_fname).read_text())


def init_dashboard():
    n = 2  # change this to control the number of HDX measurements added
    input_control = ctrl.control_panels["PeptideRFUFileInputControl"]
    for i, (k, v) in enumerate(state_spec.items()):
        if i == n:
            break
        load_state_rfu(ctrl, v, data_dir=data_dir, name=k)
        input_control._add_single_dataset_spec()

    input_control._action_load_datasets()

    # if n > 1:
    #     diff = ctrl.control_panels['DifferentialControl']
    #     diff._action_add_comparison()


# pn.state.onload(reload_dashboard)
# pn.state.onload(reload_tables)
pn.state.onload(init_dashboard)
# pn.state.onload(init_batch)


if __name__ == "__main__":
    Path(cfg.assets_dir).mkdir(exist_ok=True, parents=True)
    pn.serve(
        tmpl,
        show=True,
        static_dirs={"pyhdx": STATIC_DIR, "assets": str(cfg.assets_dir)},
    )

elif __name__.startswith("bokeh_app"):
    Path(cfg.assets_dir).mkdir(exist_ok=True, parents=True)
    tmpl.servable()

    # panel serve dev_gui_secB.py --show --autoreload --port 5076 --static-dirs pyhdx=C:/Users/jhsmi/pp/PyHDX/pyhdx/web/static assets=C:/Users/jhsmi/.pyhdx/assets
