from pyhdx.fileIO import read_dynamx
from pyhdx.support import fmt_export
from pyhdx import PeptideMasterTable, KineticsFitting
import os
import numpy as np
import pickle


directory = os.path.dirname(__file__)
fpath = os.path.join(directory, 'test_data', 'ecSecB_apo.csv')

data_dir = r'../data'
state = 'SecB WT apo'
control = ('Full deuteration control', 0.167)

data = read_dynamx(fpath)

pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pf.set_control(control)
states = pf.groupby_state()
series = states[state]
series.make_uniform()

temperature, pH = 273.15 + 30, 8.
kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

fr1 = kf.weighted_avg_fit()
out1 = fr1.output

fr_pfact = kf.global_fit_new(out1, use_kint=True, l1=20)
output = fr_pfact.output

fmt, hdr = fmt_export(output)
np.savetxt(os.path.join(directory, 'test_data', 'test123.txt'), output, fmt=fmt, header=hdr)
with open(os.path.join(directory, 'test_data', 'test123.pick'), 'wb') as f:
    pickle.dump(fr_pfact, f)

