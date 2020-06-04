import pyhdx
from pyhdx.support import np_from_txt
import numpy as np
import pickle

filename = r"SecB/ecSecB apo.csv"
drop_first = 1  # Number of N terminal residues to ignore
ignore_prolines = True
control_100 = ('Full deuteration control', 0.167)

series_name = 'SecB WT apo'

data = pyhdx.read_dynamx(filename)
pf = pyhdx.PeptideMasterTable(data, drop_first=drop_first, ignore_prolines=ignore_prolines)
pf.set_control(control_100)
states = pf.groupby_state()
series = states[series_name]
series.make_uniform()

with open('SecB/SecB_Apo_fd_fit1.pick', 'rb') as f:
    fit1_result = pickle.load(f)

with open('SecB/SecB_Apo_fd_fit2.pick', 'rb') as f:
    fit2_result = pickle.load(f)


