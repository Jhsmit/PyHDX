"""
Reload SecB and fitted data and launch  GUI
Run local_cluster.py in anothor thread

"""


import sys
from pathlib import Path

import pandas as pd
import panel as pn
import yaml

from pyhdx.batch_processing import StateParser
from pyhdx.fileIO import csv_to_dataframe, load_fitresult
from pyhdx.web.apps import main_app
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.utils import fix_multiindex_dtypes, load_state
from pyhdx.config import cfg, reset_config


CONFIG_NAME = "config_beta.yaml"
cfg_pth = Path().home() / ".pyhdx" / CONFIG_NAME

cfg.load_config(cfg_pth)

# def printfunc(args):
# 1/0

# sys.stdout = printfunc


sys._excepthook = sys.excepthook


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


ctrl, tmpl = main_app()

cwd = Path(__file__).parent
root_dir = cwd.parent.parent

web_data_dir = root_dir / "tests" / "test_data" / "output" / "web"
input_data_dir = root_dir / "tests" / "test_data" / "input"


filenames = ["ecSecB_apo.csv", "ecSecB_dimer.csv"]
file_dict = {fname: (input_data_dir / fname).read_bytes() for fname in filenames}

# batch_fname = 'data_states_deltas.yaml' # secb apo / dimer but artificial delta C/N tail, needs additional files
batch_fname = "data_states.yaml"  #


hdx_spec = yaml.safe_load(Path(input_data_dir / batch_fname).read_text())
# pdb_string = (web_data_dir / '1qyn.pdb').read_text()
pdb_string = (web_data_dir / "5JTR_mod.pdb").read_text()


def reload_tables():
    src = ctrl.sources["main"]
    src.add_table("peptides", csv_to_dataframe(web_data_dir / "peptides.csv"))
    src.add_table("d_uptake", csv_to_dataframe(web_data_dir / "d_uptake.csv"))
    # src.add_table("rfu", csv_to_dataframe(web_data_dir / "rfu.csv"))  # is missing rfu sd's
    src.add_table("dG", csv_to_dataframe(web_data_dir / "dG.csv"))
    src.add_table("ddG_comparison", csv_to_dataframe(web_data_dir / "ddG_comparison.csv"))
    src.add_table("rates", csv_to_dataframe(web_data_dir / "rates.csv"))
    src.param.trigger("updated")

    # ctrl.views['protein'].object = pdb_string


def reload_dashboard():
    source = ctrl.sources["dataframe"]
    for ds in [
        "peptides",
        "peptides_mse",
        "d_calc",
        "rfu",
        "rates",
        "global_fit",
        "losses",
    ]:
        df = csv_to_dataframe(web_data_dir / f"{ds}.csv")
        source.add_df(df, ds)

    # Temporary workaround for comment characters in csv files
    ds = "colors"
    df = pd.read_csv(web_data_dir / f"{ds}.csv", header=[0, 1, 2], index_col=0, skiprows=3)
    source.add_df(df, ds)


def init_batch():
    input_control = ctrl.control_panels["PeptideFileInputControl"]
    input_control.input_mode = "Batch"

    file_dict = {fname: (input_data_dir / fname).read_bytes() for fname in filenames}
    input_control.widgets["input_files"].filename = list(file_dict.keys())
    input_control.input_files = list(file_dict.values())

    input_control.batch_file = Path(input_data_dir / batch_fname).read_bytes()
    input_control._action_load_datasets()

    fit_control = ctrl.control_panels["FitControl"]

    fit_control.r1 = 0.05
    fit_control.r2 = 0.1
    fit_control.epochs = 200
    fit_control.stop_loss = 0.001
    fit_control.patience = 10000
    fit_control.learning_rate = 100


def init_dashboard():
    n = 2  # change this to control the number of HDX measurements added
    input_control = ctrl.control_panels["PeptideFileInputControl"]
    # states = ["SecB_tetramer"]
    states = ["SecB_tetramer", "SecB_dimer"]
    load_state(input_control, hdx_spec, data_dir=input_data_dir, states=states)

    input_control._action_load_datasets()

    d_uptake_control = ctrl.control_panels["DUptakeFitControl"]
    d_uptake_control.repeats = 2

    ctrl.sources["pdb"].add_from_string(pdb_string, "1qyn")

    src = ctrl.sources["main"]
    # df = csv_to_dataframe(web_data_dir / 'd_uptake.csv')
    # df.columns = fix_multiindex_dtypes(df.columns)
    # src.add_table('d_uptake', df)
    # src.updated = True
    #
    # df = csv_to_dataframe(web_data_dir / 'dG.csv')
    # df.columns = fix_multiindex_dtypes(df.columns)
    # src.add_table('dG', df)

    # src.updated = True

    # guess_control = ctrl.control_panels['InitialGuessControl']
    # guess_control._action_fit()
    #
    # fit_control = ctrl.control_panels['FitControl']
    #
    # fit_control.r1 = 0.05
    # fit_control.r2 = 0.1
    # fit_control.epochs = 200
    # fit_control.stop_loss = 0.001
    # fit_control.patience = 10000
    # fit_control.learning_rate = 100
    #
    # pdbe = ctrl.views['protein']
    # #ctrl.views['protein'].object = pdb_string
    # #
    # fit_result = load_fitresult(fitresult_dir)
    # src = ctrl.sources['main']
    # src.add(fit_result, 'fit_1')
    #
    # if n > 1:
    #     diff = ctrl.control_panels['DifferentialControl']
    #     diff._action_add_comparison()

    # if n > 1:
    #     diff = ctrl.control_panels['DifferentialControl']
    #     diff._action_add_comparison()


def init_manual():
    input_control = ctrl.control_panels["PeptideFileInputControl"]

    file_dict = {fname: (input_data_dir / fname).read_bytes() for fname in filenames}
    input_control.widgets["input_files"].filename = list(file_dict.keys())
    input_control.input_files = list(file_dict.values())

    input_control.fd_state = "Full deuteration control"
    input_control.fd_exposure = 10.020000000000001

    input_control.exp_state = "SecB WT apo"


# pn.state.onload(reload_dashboard)
# pn.state.onload(reload_tables)
pn.state.onload(init_dashboard)
# pn.state.onload(init_batch)
# pn.state.onload(init_manual)


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
