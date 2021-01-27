from pathlib import Path
import numpy as np
from pyhdx import PeptideMasterTable, KineticsFitting, read_dynamx
from pyhdx.fileIO import txt_to_np, fmt_export
import pickle

fit_dir = 'test_data'
directory = os.path.dirname(__file__)
np.random.seed(43)

fpath = os.path.join(directory, 'test_data', 'ecSecB_apo.csv')
data = read_dynamx(fpath)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(('Full deuteration control', 0.167))
states = pmt.groupby_state()

series = states['SecB WT apo']


