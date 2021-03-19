"""
Reload SecB and fitted data and launch  GUI

"""

from pyhdx.fileIO import read_dynamx, txt_to_np, csv_to_dataframe
from pyhdx import PeptideMasterTable
import pickle
from pyhdx.panel.apps import _main_app
from pyhdx.panel.utils import reload_previous
from pyhdx.panel.base import DEFAULT_COLORS, STATIC_DIR
from pyhdx.panel.data_sources import DataSource
import panel as pn
import numpy as np
from pathlib import Path

tmpl, ctrl = _main_app()
directory = Path(__file__).parent
fpath = directory / 'test_data' / 'ecSecB_apo.csv'

dic = {}
dic['file_path'] = directory / 'test_data' / 'ecSecB_apo.csv'
dic['norm_mode'] = 'Exp'
dic['fd_state'] = 'Full deuteration control'
dic['fd_exposure'] = 0.167
dic['exp_state'] = 'SecB WT apo'

src_file = directory / 'test_data' / 'ecSecB_torch_fit.txt'
df = csv_to_dataframe(src_file)
data_dict = df.to_dict(orient='series')

print(data_dict)


data_dict['color'] = np.full_like(data_dict['r_number'], fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')
data_source = DataSource(data_dict, x='r_number', tags=['mapping', 'pfact', 'deltaG'],
                         renderer='circle', size=10, name='global_fit')

dic['sources'] = {}
#dic['sources']['global_fit'] = data_source  #todo: on_load!

dic['rcsb_id'] = '1qyn'

cluster = None
ctrl = reload_previous(dic, ctrl)


if __name__ == '__main__':
    pn.serve(tmpl, show=False, static_dirs={'pyhdx': STATIC_DIR})
