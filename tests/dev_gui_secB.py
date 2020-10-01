"""
Reload SecB and fitted data and launch  GUI

"""


from pyhdx.fileIO import read_dynamx
from pyhdx import PeptideMasterTable
from pyhdx.support import np_from_txt
import pickle
import os
from pyhdx.panel.apps import _main_app
from pyhdx.panel.utils import reload_previous
from pyhdx.panel.base import DEFAULT_COLORS
from pyhdx.panel.data_sources import DataSource
import panel as pn
import numpy as np


tmpl, ctrl = _main_app()
directory = os.path.dirname(__file__)

dic = {}
dic['file_path'] = os.path.join(directory, 'test_data', 'ecSecB_apo.csv')
dic['norm_mode'] = 'Exp'
dic['norm_state'] = 'Full deuteration control'
dic['norm_exposure'] = 0.167
dic['exp_state'] = 'SecB WT apo'

src_file = os.path.join(directory, 'test_data', 'SecB WT apo_pfact_linear.txt')
array = np_from_txt(src_file)
data_dict = {name: array[name] for name in array.dtype.names}


data_dict['color'] = np.full_like(array, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')
data_dict['color'][np.isnan(data_dict['log_P'])] = np.nan
data_dict['pfact'] = 10**data_dict['log_P']

data_source = DataSource(data_dict, x='r_number', y='pfact', tags=['mapping', 'pfact'],
                         renderer='circle', size=10)

dic['sources'] = {}
dic['sources']['pfact'] = data_source

dic['rcsb_id'] = '1qyn'

cluster = None
ctrl = reload_previous(dic, ctrl)


if __name__ == '__main__':
    pn.serve(tmpl, show=False)

