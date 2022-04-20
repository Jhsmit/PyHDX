"""
Reload folding data and GUI
Run local_cluster.py in anothor thread

"""

import sys
from pathlib import Path

import pandas as pd
import panel as pn
import yaml

from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.fileIO import csv_to_dataframe, load_fitresult
from pyhdx.fileIO import csv_to_protein
from pyhdx.web.apps import rfu_app
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.utils import load_state, fix_multiindex_dtypes

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


cwd = Path(__file__).parent.resolve()
root_dir = cwd.parent.parent
test_dir = cwd / "test_data_folding"

# fpath_1 = root_dir / 'tests' / 'test_data' / 'ecSecB_apo.csv'
# fpath_2 = root_dir / 'tests' / 'test_data' / 'ecSecB_dimer.csv'
# fitresult_dir = root_dir / 'tests' / 'test_data' / 'output' / 'ecsecb_tetramer_dimer'
#
#
# data_dir = root_dir / 'tests' / 'test_data' / 'input'
# #data_dir = cwd / 'rinky_data'
#
# yaml_dict = yaml.safe_load(Path(data_dir / 'data_states.yaml').read_text())
# pdb_string = (test_dir / '1qyn.pdb').read_text()


def reload_tables():
    test_dir = cwd / "test_data"
    src = ctrl.sources["main"]

    df = csv_to_dataframe(test_dir / "peptides.csv")
    # names = df.columns.names
    # df = df.convert_dtypes()
    # df.columns.names = names
    table_names = [
        "rfu_residues.csv",
        "rates.csv",
        "peptides.csv",
        "dG_fits.csv",
        "ddG_comparison.csv",
        "d_calc.csv",
        "loss.csv",
        "peptide_mse.csv",
    ]
    for name in table_names:
        try:
            df = csv_to_dataframe(test_dir / name)
            df.columns = fix_multiindex_dtypes(df.columns)
            src.tables[name.split(".")[0]] = df
        except Exception as e:
            print(e)
            print("not loaded:", name)

    src.param.trigger("updated")

    # ctrl.views['protein'].object = pdb_string


def init_mbp():
    file_input = ctrl.control_panels["PeptideRFUFileInputControl"]

    # -------------------------------------- #
    filename = "MBPwt_4C_folding.csv"
    binary_data = (test_dir / filename).read_bytes()
    file_input.input_files = [binary_data]

    file_input.fd_state = "MBP_wt_native"
    file_input.fd_exposure = 960.000061 * 60

    file_input.nd_state = "FD"
    file_input.nd_exposure = 0.001000 * 60

    file_input.exp_state = "MBP_4C_fold_kinetix"
    file_input.exp_exposures = file_input.exp_exposures[1:]

    file_input.c_term = 375

    file_input._action_add_dataset()

    # -------------------------------------- #
    filename = "MBPp101_4C_fold.csv"
    binary_data = (test_dir / filename).read_bytes()
    file_input.input_files = [binary_data]

    file_input.fd_state = "NATIVE_MBPp101"
    file_input.fd_exposure = 60.000004 * 60

    file_input.nd_state = "FD_MBPp101"
    file_input.nd_exposure = 0.001000 * 60

    file_input.exp_state = "MBPp101_4C_folding_kinetix"
    file_input.exp_exposures = file_input.exp_exposures[1:]

    file_input.c_term = 375

    file_input._action_add_dataset()

    diff = ctrl.control_panels["DifferentialControl"]
    diff._action_add_comparison()

    pdb_src = ctrl.sources["pdb"]
    pdb_src.add_from_pdb("1MPD")


def init_ppia():
    filename = "wt_ppiA_folding_4Cmodif_230120.csv"
    binary_data = (test_dir / filename).read_bytes()

    file_input = ctrl.control_panels["PeptideRFUFileInputControl"]
    file_input.input_files = [binary_data]

    file_input.fd_state = "Native folded"
    file_input.fd_exposure = 3600.00024

    file_input.nd_state = "FD"
    file_input.nd_exposure = 0.06

    file_input.exp_state = "folding_4C_10secLabelling"
    file_input.exp_exposures = file_input.exp_exposures[1:]

    # file_input._action_add_dataset()

    pdb_src = ctrl.sources["pdb"]
    pdb_src.add_from_pdb("1qyn")
    # pdb_src.add_from_string(pdb_#string, '1qyn')


# if __name__ == '__main__':
# pn.state.onload(reload_dashboard)
# pn.state.onload(reload_tables)
pn.state.onload(init_mbp)

if __name__ == "__main__":
    pn.serve(tmpl, show=True, static_dirs={"pyhdx": STATIC_DIR})

elif __name__.startswith("bokeh_app"):
    tmpl.servable()

# ctrl.template.servable()


# panel serve --show --autoreload --static-dirs pyhdx=C:\Users\jhsmi\pp\PyHDX\pyhdx\web\static
