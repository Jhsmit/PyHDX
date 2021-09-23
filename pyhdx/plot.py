"""
Outdated module
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import proplot as pplt
import pyhdx
from pyhdx.support import autowrap, rgb_to_hex
from pyhdx.fileIO import load_fitresult
import warnings


no_coverage = '#8c8c8c'
node_pos = [10, 25, 40]  # in kJ/mol
linear_colors = ['#ff0000', '#00ff00', '#0000ff']  # red, green, blue
rgb_norm = plt.Normalize(node_pos[0], node_pos[-1], clip=True)
rgb_cmap = mpl.colors.LinearSegmentedColormap.from_list("rgb_cmap", list(zip(rgb_norm(node_pos), linear_colors)))
rgb_cmap.set_bad(color=no_coverage)

diff_colors = ['#54278e', '#ffffff', '#006d2c'][::-1]
diff_node_pos = [-10, 0, 10]
diff_norm = plt.Normalize(diff_node_pos[0], diff_node_pos[-1], clip=True)
diff_cmap = mpl.colors.LinearSegmentedColormap.from_list("diff_cmap", list(zip(diff_norm(diff_node_pos), diff_colors)))
diff_cmap.set_bad(color=no_coverage)

cbar_width = 0.075

dG_ylabel = 'ΔG (kJ/mol)'
ddG_ylabel = 'ΔΔG (kJ/mol)'

r_xlabel = 'Residue Number'


errorbar_kwargs = {
    'fmt': 'o',
    'ecolor': 'k',
    'elinewidth': 0.3,
    'markersize': 0,
    'alpha': 0.75
}

scatter_kwargs = {
    's': 7
}

def plot_residue_map(pm, scores=None, ax=None, cmap='jet', bad='k', cbar=True, **kwargs): # pragma: no cover
    """
    FUNCTION IS MOST LIKELY OUT OF DATE

    Parameters
    ----------
    pm
    scores
    ax
    cmap
    bad
    cbar
    kwargs

    Returns
    -------

    """

    warnings.warn("This function will be removed", DeprecationWarning)

    img = (pm.X > 0).astype(float)
    if scores is not None:
        img *= scores[:, np.newaxis]
    elif pm.rfu is not None:
        img *= pm.rfu[:, np.newaxis]

    ma = np.ma.masked_where(img == 0, img)
    cmap = mpl.cm.get_cmap(cmap)
    cmap.set_bad(color=bad)

    ax = plt.gca() if ax is None else ax
    ax.set_facecolor(bad)

    im = ax.imshow(ma, cmap=cmap, **kwargs)
    if cbar:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Uptake (%)')

    ax.set_xlabel('Residue number')
    ax.set_ylabel('Peptide index')







def plot_peptides(pm, ax, wrap=None,
                  color=True, labels=False, cbar=False,
                  intervals='corrected', cmap='jet', **kwargs):
    """

    TODO: needs to be checked if intervals (start, end) are still accurately taking inclusive, exclusive into account
    Plots peptides as rectangles in the provided axes

    Parameters
    ----------
    pm
    wrap
    ax
    color
    labels
    cmap
    kwargs

    Returns
    -------

    """

    wrap = wrap or autowrap(pm.data['start'], pm.data['end'])
    rect_kwargs = {'linewidth': 1, 'linestyle': '-', 'edgecolor': 'k'}
    rect_kwargs.update(kwargs)

    cmap = mpl.cm.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    i = -1

    for p_num, idx in enumerate(pm.data.index):
        e = pm.data.loc[idx]
        if i < -wrap:
            i = -1

        if color:
            c = cmap(norm(e['rfu']))
        else:
            c = '#707070'

        if intervals == 'corrected':
            start, end = 'start', 'end'
        elif intervals == 'original':
            start, end = '_start', '_end'
        else:
            raise ValueError(f"Invalid value '{intervals}' for keyword 'intervals', options are 'corrected' or 'original'")

        width = e[end] - e[start]
        rect = Rectangle((e[start] - 0.5, i), width, 1, facecolor=c, **rect_kwargs)
        ax.add_patch(rect)
        if labels:
            rx, ry = rect.get_xy()
            cy = ry
            cx = rx
            ax.annotate(str(p_num), (cx, cy), color='k', fontsize=6, va='bottom', ha='right')

        i -= 1

    if cbar:
        scalar_mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(scalar_mappable, label='Percentage D')

    ax.set_ylim(-wrap, 0)
    end = pm.interval[1]
    ax.set_xlim(0, end)
    ax.set_yticks([])


def plot_fitresults(fitresult_path, plots='all', renew=False):
    #fit_result = csv_to_dataframe(fitresult_path / 'fit_result.csv')

    history_path = fitresult_path / 'model_history.csv'
    check_exists = lambda x: False if renew else x.exists()
    try: # temp hack as batch results do not store hdxms
        fit_result = load_fitresult(fitresult_path)
        df = fit_result.output

        dfs = [df]
        names = ['']
        hdxm_s = [fit_result.data_obj]
        loss_list = [fit_result.losses]
        if history_path.exists():
            history_list = [csv_to_dataframe(history_path)]
        else:
            history_list = []

    except FileNotFoundError:
        df = csv_to_dataframe(fitresult_path / 'fit_result.csv')
        dfs = [df[c] for c in df.columns.levels[0]]
        names = [c + '_' for c in df.columns.levels[0]]
        loss_list = [csv_to_dataframe(fitresult_path / 'losses.csv')]

        hdxm_s = []

        if history_path.exists():
            history_df = csv_to_dataframe(history_path)
            history_list = [history_df[c] for c in history_df.columns.levels[0]]
        else:
            history_list = []

    full_width = 170 / 25.4
    width = 120 / 25.4
    aspect = 4
    cmap = rgb_cmap
    norm = rgb_norm

    COV_SCALE = 1.

    if plots == 'all':
        plots = ['losses', 'deltaG', 'pdf', 'coverage', 'history']

    if 'losses' in plots:
        for loss_df in loss_list:  # Mock loop to use break
            output_path = fitresult_path / 'losses.png'
            if check_exists(output_path):
                break

#            losses = loss_df.drop('reg_percentage', axis=1)
            loss_df.plot()

            mse_loss = loss_df['mse_loss']
            reg_loss = loss_df.iloc[:, 1:].sum(axis=1)
            reg_percentage = 100*reg_loss / (mse_loss + reg_loss)
            fig = plt.gcf()
            ax = plt.gca()
            ax1 = ax.twinx()
            reg_percentage.plot(ax=ax1, color='k')
            ax1.set_xlim(0, None)
            plt.savefig(output_path)
            plt.close(fig)

    if 'deltaG' in plots:
        for result, name in zip(dfs, names):
            output_path = fitresult_path / f'{name}deltaG.png'
            if check_exists(output_path):
                break

            fig, axes = pplt.subplots(nrows=1, width=width, aspect=aspect)
            ax = axes[0]

            yvals = result['deltaG'] * 1e-3
            rgba_colors = cmap(norm(yvals), bytes=True)
            hex_colors = rgb_to_hex(rgba_colors)
            ax.scatter(result.index, yvals, c=hex_colors, **scatter_kwargs)
            ylim = ax.get_ylim()
            ax.errorbar(result.index, yvals, yerr=result['covariance'] * 1e-3 * COV_SCALE, **errorbar_kwargs, zorder=-1)

            ax.format(ylim=ylim, ylabel=dG_ylabel, xlabel=r_xlabel)

            plt.savefig(output_path, transparent=False)
            plt.close(fig)

    if 'pdf' in plots:
        for i in range(1):
            output_path = fitresult_path / 'fit_report'
            if check_exists(fitresult_path / 'fit_report.pdf'):
                break

            output = pyhdx.Output(fit_result)

            report = pyhdx.Report(output, title=f'Fit report {fit_result.data_obj.name}')
            report.add_peptide_figures()
            report.generate_pdf(output_path)

    if 'coverage' in plots:
        for hdxm in hdxm_s:
            output_path = fitresult_path / f'{hdxm.name}_coverage.png'
            if check_exists(output_path):
                break

            n_rows = int(np.ceil(len(hdxm.timepoints) / 2))

            fig, axes = pplt.subplots(ncols=2, nrows=n_rows, sharex=True, width=full_width, aspect=4)
            axes_list = list(axes[:, 0]) + list(axes[:, 1])

            for label, ax, pm in zip(hdxm.timepoints, axes_list, hdxm):
                plot_peptides(pm, ax, linewidth=0.5)
                ax.format(title=label, xlabel=r_xlabel)

            plt.savefig(output_path, transparent=False)
            plt.close(fig)

    if 'history' in plots:
        for h_df, name in zip(history_list, names):
            output_path = fitresult_path / f'{name}history.png'
            if check_exists(output_path):
                break

            num = len(h_df.columns)
            max_epochs = max([int(c) for c in h_df.columns])

            cmap = mpl.cm.get_cmap('winter')
            norm = mpl.colors.Normalize(vmin=1, vmax=max_epochs)
            colors = iter(cmap(np.linspace(0, 1, num=num)))

            fig, axes = pplt.subplots(nrows=1, width=width, aspect=aspect)
            ax = axes[0]
            for key in h_df:
                c = next(colors)
                to_hex(c)

                ax.scatter(h_df.index, h_df[key] * 1e-3, color=to_hex(c), **scatter_kwargs)
            ax.format(xlabel=r_xlabel, ylabel=dG_ylabel)

            values = np.linspace(0, max_epochs, endpoint=True, num=num)
            colors = cmap(norm(values))
            tick_labels = np.linspace(0, max_epochs, num=5)

            cbar = fig.colorbar(colors, values=values, ticks=tick_labels, space=0, width=cbar_width, label='Epochs')
            ax.format(yticklabelloc='None', ytickloc='None')

            plt.savefig(output_path)
            plt.close(fig)