"""
Reload SecB and fitted data and launch  GUI

"""

from pyhdx.fileIO import read_dynamx, txt_to_np, csv_to_dataframe
from pyhdx import PeptideMasterTable
import pickle
from pyhdx.panel.apps import main_app
from pyhdx.panel.utils import reload_previous
from pyhdx.panel.base import DEFAULT_COLORS, STATIC_DIR
from pyhdx.panel.sources import DataSource
import panel as pn
import numpy as np
from pathlib import Path

#temporary imports
from pyhdx.support import rgb_to_hex
import matplotlib.pyplot as plt


ctrl = main_app()
directory = Path(__file__).parent
fpath_1 = directory / 'test_data' / 'ecSecB_apo.csv'
fpath_2 = directory / 'test_data' / 'ecSecB_dimer.csv'

fpaths = [fpath_1, fpath_2]
files = [p.read_bytes() for p in fpaths]

file_input = ctrl.control_panels['PeptideFileInputControl']

file_input.input_files = files
file_input.fd_state = 'Full deuteration control'
file_input.fd_exposure = 0.167

file_input.exp_state = 'SecB WT apo'
file_input.dataset_name = 'testname_123'
file_input._action_add_dataset()

file_input.exp_state = 'SecB his dimer apo'
file_input.dataset_name = 'SecB his dimer apo'  # todo catch error duplicate name
file_input._action_add_dataset()

initial_guess = ctrl.control_panels['InitialGuessControl']
initial_guess._action_fit()

fit_control = ctrl.control_panels['FitControl']
fit_control.epochs = 10

fit_control._do_fitting()

classification = ctrl.control_panels['ClassificationControl']
classification.widgets['select_1'].value = '*'
classification.widgets['select_2'].value = 'deltaG'

classification.mode = 'Continuous'
classification._action_linear()
classification.color_set_name = 'colorset test'
classification._action_add_colorset()


file_export = ctrl.control_panels['FileExportControl']
sio = file_export.table_export_callback()



#
if __name__ == '__main__':
    pn.serve(ctrl.template, show=True
             , static_dirs={'pyhdx': STATIC_DIR})
