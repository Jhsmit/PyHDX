import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def plot_residue_map(pm, scores=None, ax=None, cmap='jet', bad='k', cbar=True, **kwargs):
    img = (pm.big_X > 0).astype(float)
    if scores is not None:
        img *= scores[:, np.newaxis]
    elif pm.scores is not None:
        img *= pm.scores[:, np.newaxis]

    ma = np.ma.masked_where(img == 0, img)
    cmap = matplotlib.cm.get_cmap(cmap)
    cmap.set_bad(color=bad)

    ax = plt.gca() if ax is None else ax
    ax.set_facecolor(bad)

    im = ax.imshow(ma, cmap=cmap, **kwargs)
    if cbar:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Uptake (%)')

    ax.set_xlabel('Residue number')
    ax.set_ylabel('Peptide index')
