"""Load a fit result from a a directory where a fit result was saved with save_fitresult"""
from pyhdx.fileIO import load_fitresult
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


current_dir = Path().cwd()
fit_result = load_fitresult(current_dir / 'output' / 'SecB_fit')


time = np.logspace(-3, 2, num=100)

d_calc = fit_result(time)
d_exp = fit_result.hdxm_set.d_exp

i = 20  # index of the protein to view

fit_result.losses[['total_loss', 'mse_loss', 'reg_loss']].plot()

fig, ax = plt.subplots()
ax.scatter(fit_result.hdxm_set.timepoints, d_exp[i], color='k')
ax.plot(time, d_calc[i], color='r')
ax.set_xscale('log')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Corrected D-uptake (Da)')
ax.set_title(f'Peptide index: {i}')

fig, ax = plt.subplots()
ax.set_aspect(2)
ax.scatter(fit_result.output.index, fit_result.output['deltaG']*1e-3)
ax.set_xlabel('Residue number')
ax.set_ylabel('ΔG (kJ/mol)')
plt.show()