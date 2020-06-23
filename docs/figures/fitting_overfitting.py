from pyhdx import read_dynamx, PeptideMasterTable, KineticsFitting
from pyhdx.support import get_reduced_blocks, get_original_blocks
from pyhdx.plot import plot_peptides
import matplotlib.pyplot as plt

import os
import numpy as np

data_dir = '../../tests/test_data'
filename = 'ecSecB_apo.csv'
refit = False

fpath = os.path.join(data_dir, filename)

data = read_dynamx(fpath)
master_table = PeptideMasterTable(data, drop_first=0, ignore_prolines=False)
master_table.set_control(('Full deuteration control', 0.167))

states = master_table.groupby_state()
print(states.keys())
series = states['SecB WT apo']
series.make_uniform()
split = series.split()
key = list(split)[1]
series = split[key]

kf = KineticsFitting(series, bounds=(0, 200))
print(kf.bounds)





if refit:
    fr1 = kf.weighted_avg_fit()
    arr = fr1.get_output(['tau1', 'tau2', 'r'])

    fr_orig = kf.blocks_fit(arr, block_func=get_original_blocks)
    out_orig = fr_orig.get_output(['rate', 'tau1', 'tau2'])
    np.save('fr_orig.npy', out_orig)

    fr_red = kf.blocks_fit(arr)
    out_red = fr_red.get_output(['rate', 'tau1', 'tau2'])
    np.save('fr_red.npy', out_red)

else:
    out_orig = np.load('fr_orig.npy')
    out_red = np.load('fr_red.npy')


fig, ax = plt.subplots(1, figsize=(12, 3.5))
ax.set_yscale('log')
ax.scatter(out_orig['r_number'], out_orig['rate'], color='r', label='original', marker='v')
ax.scatter(out_red['r_number'], out_red['rate'], color='b', label='reduced')
ax.legend()
ax.set_xlabel('Residue')
ax.set_ylabel('Rate (min⁻¹)')
ax.set_xlim(16, 59)

cov = series.cov

plt.tight_layout()
#plt.show()
plt.savefig('overfitting.png')


