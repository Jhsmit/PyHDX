"""
Reload SecB and fitted data and launch  GUI

"""

from pyhdx.fileIO import read_dynamx, txt_to_np
from pyhdx import PeptideMasterTable
import pickle
from pyhdx.panel.apps import _color_matrix_app
from pyhdx.panel.utils import reload_previous
from pyhdx.panel.base import DEFAULT_COLORS, STATIC_DIR
from pyhdx.panel.sources import DataSource
import panel as pn
import numpy as np
from pathlib import Path

tmpl, ctrl = _color_matrix_app()
directory = Path(__file__).parent

input_file_path = directory / 'test_data' / 'ec_SecB_scores.txt'

file_input_ctrl = ctrl.control_panels['MatrixMappingFileInputControl']
file_input_widget = file_input_ctrl.widget_dict['input_file']
file_input_widget.filename = input_file_path.name
file_input_widget.value = input_file_path.read_bytes()

#
# fpath = directory / 'test_data' / 'ecSecB_apo.csv'
#
# dic = {}
# dic['file_path'] = directory / 'test_data' / 'ecSecB_apo.csv'
# dic['norm_mode'] = 'Exp'
# dic['fd_state'] = 'Full deuteration control'
# dic['fd_exposure'] = 0.167
# dic['exp_state'] = 'SecB WT apo'
#
# file_input = ctrl.control_panels['PeptideFileInputControl']
#
#
#
# src_file = directory / 'test_data' / 'ecSecB_torch_fit.txt'
# array = txt_to_np(src_file)
# data_dict = {name: array[name] for name in array.dtype.names}

#
# data_dict['color'] = np.full_like(array, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')
# data_source = DataSource(data_dict, x='r_number', tags=['mapping', 'pfact', 'deltaG'],
#                          renderer='circle', size=10, name='global_fit')
#
# dic['sources'] = {}
# #dic['sources']['global_fit'] = data_source  #todo: on_load!
#
# dic['rcsb_id'] = '1qyn'
#
# cluster = None
# ctrl = reload_previous(dic, ctrl)
#
#
if __name__ == '__main__':
    pn.serve(tmpl, show=False, static_dirs={'pyhdx': STATIC_DIR})
