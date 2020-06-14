from pyhdx import read_dynamx, PeptideMasterTable
from pyhdx.support import get_reduced_blocks
from pyhdx.plot import plot_peptides
import matplotlib.pyplot as plt

import os
import numpy as np

data_dir = '../../tests/test_data'
filename = 'ecSecB_apo.csv'
fpath = os.path.join(data_dir, filename)

data = read_dynamx(fpath)
master_table = PeptideMasterTable(data, drop_first=0, ignore_prolines=False)

states = master_table.groupby_state()
print(states.keys())
series = states['SecB WT apo']
split = series.split()
key = list(split)[1]
cov = split[key].cov


def add_blocks(ax, positions, color):
    for pos in positions:
        ax.plot([pos, pos], [-40, 2], color=color, linewidth=2)  # linestyle=(0, (1, 1))

    text_x = positions[:-1] + np.diff(positions) / 2
    for i, x in enumerate(text_x[text_x < 58]):
        ax.text(x, -30, f'$b_{{{i + 1}}}$', horizontalalignment='center', fontsize=12)

    for i in range(cov.end + 5):
        ax.axvline(i + 0.5, alpha=0.2, zorder=-3, color='grey')

    ax.set_ylabel('Peptides')
    ax.set_xlabel('Residue number')
    ax.set_xlim(16, 59)
    ax.set_ylim(-32, 2)

    plt.tight_layout()


fig, ax = plt.subplots(1, figsize=(12, 3.5))
plot_peptides(cov, len(cov), ax, color=False)

positions = np.cumsum(cov.block_length) + cov.start
positions = np.insert(positions, 0, cov.start) - 0.5
add_blocks(ax, positions, 'r')
plt.savefig('Blocks_original.png', dpi=300)

fig, ax = plt.subplots(1, figsize=(12, 3.5))
plot_peptides(cov, len(cov), ax, color=False)
red_blocks = get_reduced_blocks(cov)

positions = np.cumsum(red_blocks) + cov.start
positions = np.insert(positions, 0, cov.start) - 0.5
rect = ax.patches[6].set_color('r')
add_blocks(ax, positions, 'b')
plt.savefig('Blocks_reduced', dpi=300)
