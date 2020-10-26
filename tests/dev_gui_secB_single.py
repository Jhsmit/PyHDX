"""
Reload SecB and fitted data and launch single GUI

"""
from pyhdx.fileIO import read_dynamx, txt_to_np
from pyhdx import PeptideMasterTable
import pickle
import os
from pyhdx.panel.apps import _single_app
from pyhdx.panel.utils import reload_previous
from pyhdx.panel.base import DEFAULT_COLORS
from pyhdx.panel.data_sources import DataSource
import panel as pn
import numpy as np
from pathlib import Path

tmpl, ctrl = _single_app()
directory = Path(__file__).parent

fpath = directory / 'test_data' / 'SecB WT apo_pfact_linear.txt'
with open(fpath, 'rb') as f_obj:
    file_binary = f_obj.read()

f_input = ctrl.control_panels['SingleMappingFileInputControl']
f_input._widget_dict['input_file'].filename = str(fpath)
f_input.input_file = file_binary

f_input.dataset_name = 'DS1'
f_input._action_add_dataset()

# s_ctrl = ctrl.control_panels['SingleControl']
# s_ctrl.dataset_name = 'DS1_deltaG'
# s_ctrl.quantity = 'deltaG'
# s_ctrl._action_add_dataset()

pv_ctrl = ctrl.control_panels['ProteinViewControl']
pv_ctrl.rcsb_id = '1qyn'



if __name__ == '__main__':
    pn.serve(tmpl, show=True)

