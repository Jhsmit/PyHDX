import os
from pyhdx import PeptideMeasurements, PeptideCSVFile
import numpy as np
from pyhdx.fileIO import check_data, write_ass, write_dexp, write_expfact, write_seq
from pyhdx.exPfact.from_config import make_config
from pyhdx.exPfact.exPfact import run_config

ds1_pth = r'../../tests/test_data/ds1.csv'

csvfile = PeptideCSVFile(ds1_pth)
#p_dict = csvfile.return_by_name('PpiA-FD', 0.167)
#print(np.unique(csvfile.data['exposure']))
p_dict = csvfile.return_by_name('PpiANative', 30.000002)
state = 'PpiANative'
data_list = [pf for pf in p_dict.values() if pf.state == state]

# last time point as 100% or not?

#write_expfact(data_list, 'test')
in_dir = os.path.dirname(__file__)
out_dir = os.path.join(in_dir, 'output')
cfg = make_config('test', in_dir, out_dir)

cfg['pfact'] = os.path.join(in_dir, 'test.pfact')

run_config(cfg)
