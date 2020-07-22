from pyhdx.fileIO import read_dynamx
from pyhdx import PeptideMasterTable
import os

directory = os.path.dirname(__file__)
fpath = os.path.join(directory, 'test_data', 'ds1.csv')

data_dir = r'../data'
state = 'PpiANative'
control = ('PpiA-FD', 0.167)

data = read_dynamx(fpath)

pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pf.set_control(control)
states = pf.groupby_state()
series = states[state]
series.make_uniform()
