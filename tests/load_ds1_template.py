from pyhdx.fileIO import read_dynamx, txt_to_protein
from pyhdx import PeptideMasterTable, KineticsFitting
import os
import numpy as np
import pickle
from pathlib import Path

directory = Path(__file__).parent
fpath = directory / 'test_data' / 'ecSecB_apo.csv'

state = 'SecB WT apo'
control = ('Full deuteration control', 0.167)

data = read_dynamx(fpath)

pf = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pf.set_control(control)
states = pf.groupby_state()
series = states[state]

temperature, pH = 273.15 + 30, 8.
kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

fr1 = kf.weighted_avg_fit()
out1 = fr1.output
out1.to_file(directory / 'test_data' / 'ecSecB_guess.txt')



#
#
# fr_pfact = kf.global_fit(out1, l1=20)
# output = fr_pfact.output
#
# fmt, hdr = fmt_export(output)
# np.savetxt(os.path.join(directory, 'test_data', 'test123.txt'), output, fmt=fmt, header=hdr)
# with open(os.path.join(directory, 'test_data', 'test123.pick'), 'wb') as f:
#     pickle.dump(fr_pfact, f)
#
