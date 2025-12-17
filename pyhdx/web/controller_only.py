# %%
# %load_ext autoreload
# %autoreload 2

# %%
from pyhdx.web.controllers import PeptideFileInputControl
from pyhdx.web.main_controllers import PyHDXController
import logging
import sys

from pyhdx.web.sources import PyHDXSource


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

logger = logging.getLogger("debug_logger")
logger.setLevel(logging.DEBUG)

main = PyHDXController(
    control_panels=[(PeptideFileInputControl, {})],
    loggers={"pyhdx": logger},
    sources={"main": PyHDXSource()},
)
# %%

input_ctrl: PeptideFileInputControl = main.control_panels["PeptideFileInputControl"]


# headless testing
# input_ctrl.dataset_id


# input_ctrl.input_mode = "Database"
# input_ctrl._action_load_datasets()

input_ctrl.panel.servable()
