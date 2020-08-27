from pyhdx.fileIO import read_dynamx
from pyhdx import PeptideMasterTable
from pyhdx.support import np_from_txt
import pickle
import os
from pyhdx.panel.main import tmpl, ctrl
from pyhdx.panel.utils import reload_previous
from pyhdx.panel.base import DEFAULT_COLORS
import panel as pn
import numpy as np

directory = os.path.dirname(__file__)

dic = {}
dic['file_path'] = os.path.join(directory, 'test_data', 'ecSecB_apo.csv')
dic['norm_mode'] = 'Exp'
dic['norm_state'] = 'Full deuteration control'
dic['norm_exposure'] = 0.167
dic['exp_state'] = 'SecB WT apo'

src_file = os.path.join(directory, 'test_data', 'SecB WT apo_pfact_linear.txt')
array = np_from_txt(src_file)
src_dict = {name: array[name] for name in array.dtype.names}
src_dict['y'] = src_dict['log_P']
src_dict['color'] = np.full_like(array, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')
src_dict['color'][np.isnan(src_dict['y'])] = np.nan
dic['sources'] = {}
dic['sources']['pfact'] = src_dict

dic['rcsb_id'] = '1qyn'

cluster = None
ctrl = reload_previous(dic, ctrl)



if __name__ == '__main__':
    pn.serve(tmpl)

