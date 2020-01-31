import matplotlib as mpl
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


def make_kinetics_figure(pm_dict, cmap='cool'):
    """

    :param pm_dict: dictionary of PeptideMeasuements
    :param times: array_like of
    :param cmap: optional string indicating which colormap to use
    :return:
    """

    """returns matplotlib figure for visualization of kinetics"""

    fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 4), sharex=True, gridspec_kw={'hspace': 0})
    times = [p.exposure for p in pm_dict.values()]
    s = pm_dict[next(iter(pm_dict))]

    norm = mpl.colors.Normalize(vmin=0, vmax=np.max(times))
    normed_times = norm(times)
    r_number = np.arange(s.start, s.stop + 1)

    colors = mpl.cm.get_cmap(cmap, len(pm_dict))(normed_times)

    for c, (k, v) in zip(colors, pm_dict.items()):
        ax1.plot(r_number, v.scores_average, color=c, marker='.', linestyle='')

    ax1.set_ylabel('Score (%)')

    ax2.plot(r_number, np.ones_like(r_number), marker='.', linestyle='', color='k')

    # lc = LineCollection(self._lc_data, colors='r')
    # ax2.add_collection(lc)
    #        ax2.axhline(self.rate_max, color='r', autolim=False)
    ax2.set_yscale('log')
    ax2.set_xlabel('Residue number')
    ax2.set_ylabel('Rate constant\n (min$^{-1})$')

    fig.align_ylabels()  # todo doesnt seem to work

    fig.subplots_adjust(right=0.85)

    cbar_ax = fig.add_axes([0.87, 0.05, 0.02, 0.9])

    cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=mpl.cm.get_cmap(cmap),
                                    norm=norm,
                                    orientation='vertical')
    cb1.set_label('Time (min)')
    #   plt.tight_layout()

    return fig, (ax1, ax2, cbar_ax)
