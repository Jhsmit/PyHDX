from pathlib import Path
import numpy as np
from pyhdx import PeptideMasterTable, KineticsFitting, read_dynamx
from pyhdx.fileIO import txt_to_protein

guess = False

root_dir = Path().resolve().parent
test_data_dir = root_dir / 'tests' / 'test_data'
input_file_path = test_data_dir / 'ecSecB_apo.csv'

data = read_dynamx(test_data_dir / 'ecSecB_apo.csv', test_data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
pmt.set_control(('Full deuteration control', 0.167))
states = pmt.groupby_state()

series = states['SecB WT apo']

temperature, pH = 273.15 + 30, 8.
kf = KineticsFitting(series, bounds=(1e-2, 800), temperature=temperature, pH=pH)

if guess:
    wt_avg_result = kf.weighted_avg_fit()
    init_guess = wt_avg_result.output
else:
    init_guess = txt_to_protein(test_data_dir / 'ecSecB_guess.txt')


from time import time
t0 = time()  # 25 secondjes (defualt settings)
fr_torch = kf.global_fit(init_guess)
t1 = time()

print(t1 - t0)

import matplotlib.pyplot as plt

print(fr_torch)
loss = fr_torch.metadata['reg_loss']
print(loss[-1])

plt.plot(loss)
plt.show()