from pyhdx.panel.utils import reload_previous
from pyhdx.support import np_from_txt
import os
import numpy as np
from pyhdx.panel.base import DEFAULT_COLORS
from pyhdx.panel.main import tmpl, ctrl
import panel as pn
import pickle

dic = {}

directory = os.path.dirname(__file__)
dic['file_path'] = os.path.join(directory, 'test_data', 'simulated_data_uptake.csv')

dic['norm_mode'] = 'Theory'
dic['be_percent'] = 0.
dic['exp_state'] = 'state1'

#todo this should be moved to reload_previous function
# src_file = os.path.join(directory, 'test_data', 'fit_simulated_pfact.txt')
# array = np_from_txt(src_file)
# src_dict = {name: array[name] for name in array.dtype.names}
# src_dict['y'] = src_dict['log_P']
# src_dict['color'] = np.full_like(array, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')
#
# dic['sources'] = {}
# dic['sources']['pfact'] = src_dict
#
# with open(os.path.join(directory, 'test_data', 'fit_simulated_pfact.pick'), 'rb') as f:
#     fit_result = pickle.load(f)
#
# dic['fit_results'] = {}
# dic['fit_results']['fr_pfact'] = fit_result


ctrl = reload_previous(dic, ctrl)

if __name__ == '__main__':
    pn.serve(tmpl)