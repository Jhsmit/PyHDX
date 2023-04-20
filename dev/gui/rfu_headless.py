import sys
from logging import StreamHandler
from pathlib import Path

import panel as pn

from pyhdx.web.apps import rfu_app
import yaml

from pyhdx.web.utils import load_state

# %%

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


# %%

cwd = Path(__file__).parent
root_dir = cwd.parent.parent
data_dir = root_dir / "tests" / "test_data" / "input"

batch_fname = "PpiX_states.yaml"  # secb apo / dimer but artificial delta C/N tail
hdx_spec = yaml.safe_load(Path(data_dir / batch_fname).read_text())

# %%
# create controller / app and load data
ctrl, tmpl = rfu_app()

file_input = ctrl.control_panels["PeptideRFUFileInputControl"]
# states = ['PpiA_Folding']
# states = ['PpiB_Folding']
states = ["PpiA_Folding", "PpiB_Folding"]
# load_state_rfu(file_input, hdx_spec, data_dir = data_dir, states=states)
#
# file_input._action_load_datasets()

print(file_input.nd_control)
# %%

# len(file_input.src.hdxm_objects)


# %%

# file_input._add_single_dataset_spec()
# %%

# input_files = [data_dir / f_dict['filename'] for f_dict in hdx_spec["data_files"].values()]
# f_bytes = [f.read_bytes() for f in input_files]
# file_input.widgets["input_files"].filename = [f.name for f in input_files]
# file_input.input_files = f_bytes
# file_input._read_files()
#
# file_input.data_files
#
#

#
# #%%
# file_input.widgets["input_files"].filename
#
# #%%
# # state_spec = hdx_spec['states']
# states = states or list(state_spec.keys())
# names = names or states
