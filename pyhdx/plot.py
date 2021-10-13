from contextlib import contextmanager
from copy import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import proplot as pplt
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from scipy.stats import kde

from pyhdx.config import cfg
from pyhdx.fileIO import load_fitresult
from pyhdx.support import autowrap, color_pymol, apply_cmap

try:
    from pymol import cmd
except ModuleNotFoundError:
    cmd = None

dG_ylabel = 'ΔG (kJ/mol)'
ddG_ylabel = 'ΔΔG (kJ/mol)'
r_xlabel = 'Residue Number'


ERRORBAR_KWARGS = {
    'fmt': 'o',
    'ecolor': 'k',
    'elinewidth': 0.3,
    'markersize': 0,
    'alpha': 0.75,
    'capthick': 0.3,
    'capsize': 0.

}

SCATTER_KWARGS = {
    's': 7
}

RECT_KWARGS = {
    'linewidth': 0.5,
    'linestyle': '-',
    'edgecolor': 'k'}

CBAR_KWARGS = {
    'space': 0,
    'width': cfg.getfloat('plotting', 'cbar_width') / 25.4,
    'tickminor': True
}


def peptide_coverage_figure(data, wrap=None, cmap='turbo', norm=None, color_field='rfu', subplot_field='exposure',
                            rect_fields=('start', 'end'), rect_kwargs=None, **figure_kwargs):

    subplot_values = data[subplot_field].unique()
    sub_dfs = {value: data.query(f'`{subplot_field}` == {value}') for value in subplot_values}

    n_subplots = len(subplot_values)

    ncols = figure_kwargs.pop('ncols', min(cfg.getint('plotting', 'ncols'), n_subplots))
    nrows = figure_kwargs.pop('nrows', int(np.ceil(n_subplots / ncols)))
    figure_width = figure_kwargs.pop('width', cfg.getfloat('plotting', 'page_width')) / 25.4
    cbar_width = figure_kwargs.pop('cbar_width', cfg.getfloat('plotting', 'cbar_width')) / 25.4
    aspect = figure_kwargs.pop('aspect', cfg.getfloat('plotting', 'peptide_coverage_aspect'))

    cmap = pplt.Colormap(cmap)
    norm = norm or pplt.Norm('linear', vmin=0, vmax=1)

    start_field, end_field = rect_fields
    if wrap is None:
        wrap = max([autowrap(sub_df[start_field], sub_df[end_field]) for sub_df in sub_dfs.values()])

    fig, axes = pplt.subplots(ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, **figure_kwargs)
    rect_kwargs = rect_kwargs or {}
    axes_iter = iter(axes)
    for value, sub_df in sub_dfs.items():
        ax = next(axes_iter)
        peptide_coverage(ax, sub_df, cmap=cmap, norm=norm, color_field=color_field, wrap=wrap, cbar=False, **rect_kwargs)
        ax.format(title=f'{subplot_field}: {value}')

    for ax in axes_iter:
        ax.axis('off')

    start, end = data[start_field].min(), data[end_field].max()
    pad = 0.05*(end-start)
    axes.format(xlim=(start-pad, end+pad), xlabel=r_xlabel)

    if not cmap.monochrome:
        cbar_ax = fig.colorbar(cmap, norm, width=cbar_width)
        cbar_ax.set_label(color_field, labelpad=-0)
    else:
        cbar_ax = None

    return fig, axes, cbar_ax


def peptide_coverage(ax, data, wrap=None, cmap='turbo', norm=None, color_field='rfu', rect_fields=('start', 'end'), labels=False, cbar=True, **kwargs):
    start_field, end_field = rect_fields
    data = data.sort_values(by=[start_field, end_field])

    wrap = wrap or autowrap(data[start_field], data[end_field])
    cbar_width = kwargs.pop('cbar_width', cfg.getfloat('plotting', 'cbar_width')) / 25.4
    rect_kwargs = {**RECT_KWARGS, **kwargs}

    cmap = pplt.Colormap(cmap)
    norm = norm or pplt.Norm('linear', vmin=0, vmax=1)
    i = -1
    for p_num, idx in enumerate(data.index):
        elem = data.loc[idx]
        if i < -wrap:
            i = -1

        if color_field is None:
            color = cmap(0.5)
        else:
            color = cmap(norm(elem[color_field]))

        width =  elem[end_field] - elem[start_field]
        rect = Rectangle((elem[start_field] - 0.5, i), width, 1, facecolor=color, **rect_kwargs)
        ax.add_patch(rect)
        if labels:
            rx, ry = rect.get_xy()
            cy = ry
            cx = rx
            ax.annotate(str(p_num), (cx, cy), color='k', fontsize=6, va='bottom', ha='right')
        i -= 1

    ax.set_ylim(-wrap, 0)
    start, end = data[start_field].min(), data[end_field].max()
    pad = 0.05*(end-start)
    ax.set_xlim(start-pad, end+pad)
    ax.set_yticks([])

    if cbar and color_field:
        cbar_ax = ax.colorbar(cmap, norm=norm, width=cbar_width)
        cbar_ax.set_label(color_field, labelpad=-0)
    else:
        cbar_ax = None

    return cbar_ax


def residue_time_scatter_figure(hdxm, field='rfu', scatter_kwargs=None, **figure_kwargs):
    """per-residue per-exposure values for field  `field` by weighted averaging """

    n_subplots = hdxm.Nt
    ncols = figure_kwargs.pop('ncols', min(cfg.getint('plotting', 'ncols'), n_subplots))
    nrows = figure_kwargs.pop('nrows', int(np.ceil(n_subplots / ncols)))
    figure_width = figure_kwargs.pop('width', cfg.getfloat('plotting', 'page_width')) / 25.4
    aspect = figure_kwargs.pop('aspect', cfg.getfloat('plotting', 'residue_scatter_aspect'))

    fig, axes = pplt.subplots(ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, **figure_kwargs)
    scatter_kwargs = scatter_kwargs or {}
    axes_iter = iter(axes)
    for hdx_tp in hdxm:
        ax = next(axes_iter)
        residue_time_scatter(ax, hdx_tp, field=field, **scatter_kwargs)
        ax.format(title=f'exposure: {hdx_tp.exposure:.1f}')

    for ax in axes_iter:
        ax.axis('off')

    axes.format(xlabel=r_xlabel)
    return fig, axes


def residue_time_scatter(ax, hdx_tp, field='rfu', **kwargs):
    scatter_kwargs = {**SCATTER_KWARGS, **kwargs}
    values = hdx_tp.weighted_average(field)
    ax.scatter(values.index, values, **scatter_kwargs)


def residue_scatter_figure(hdxm_set, field='rfu', cmap='viridis', norm=None, scatter_kwargs=None,
                               **figure_kwargs):
    n_subplots = hdxm_set.Ns
    ncols = figure_kwargs.pop('ncols', min(cfg.getint('plotting', 'ncols'), n_subplots))
    nrows = figure_kwargs.pop('nrows', int(np.ceil(n_subplots / ncols)))  #todo disallow setting rows
    figure_width = figure_kwargs.pop('width', cfg.getfloat('plotting', 'page_width')) / 25.4
    cbar_width = figure_kwargs.pop('cbar_width', cfg.getfloat('plotting', 'cbar_width')) / 25.4
    aspect = figure_kwargs.pop('aspect', cfg.getfloat('plotting', 'residue_scatter_aspect'))

    cmap = pplt.Colormap(cmap)
    if norm is None:
        tps = np.unique(np.concatenate([hdxm.timepoints for hdxm in hdxm_set]))
        tps = tps[np.nonzero(tps)]
        norm = pplt.Norm('log', vmin=tps.min(), vmax=tps.max())
    else:
        tps = np.unique(np.concatenate([hdxm.timepoints for hdxm in hdxm_set]))

    fig, axes = pplt.subplots(ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, **figure_kwargs)
    axes_iter = iter(axes)
    scatter_kwargs = scatter_kwargs or {}
    for hdxm in hdxm_set:
        ax = next(axes_iter)
        residue_scatter(ax, hdxm, cmap=cmap, norm=norm, field=field, cbar=False, **scatter_kwargs)

    for ax in axes_iter:
        ax.axis('off')

    #todo function for this?
    locator = pplt.Locator(norm(tps))
    cbar_ax = fig.colorbar(cmap, width=cbar_width, ticks=locator)
    formatter = pplt.Formatter('simple', precision=1)
    cbar_ax.ax.set_yticklabels([formatter(t) for t in tps])
    cbar_ax.set_label('Exposure time (s)', labelpad=-0)

    axes.format(xlabel=r_xlabel)

    return fig, axes, cbar_ax


def residue_scatter(ax, hdxm, field='rfu', cmap='viridis', norm=None, cbar=True, **kwargs):
    cmap = pplt.Colormap(cmap)
    tps = hdxm.timepoints[np.nonzero(hdxm.timepoints)]
    norm = norm or pplt.Norm('log', tps.min(), tps.max())

    cbar_width = kwargs.pop('cbar_width', cfg.getfloat('plotting', 'cbar_width')) / 25.4
    scatter_kwargs = {**SCATTER_KWARGS, **kwargs}
    for hdx_tp in hdxm:
        if isinstance(norm, mpl.colors.LogNorm) and hdx_tp.exposure == 0.:
            continue
        values = hdx_tp.weighted_average(field)
        color = cmap(norm(hdx_tp.exposure))
        scatter_kwargs['color'] = color
        ax.scatter(values.index, values, **scatter_kwargs)

    if cbar:
        locator = pplt.Locator(norm(tps))
        cbar_ax = ax.colorbar(cmap, width=cbar_width, ticks=locator)
        formatter = pplt.Formatter('simple', precision=1)
        cbar_ax.ax.set_yticklabels([formatter(t) for t in tps])
        cbar_ax.set_label('Exposure time (s)', labelpad=-0)


def dG_scatter_figure(data, norm=None, cmap=None, scatter_kwargs=None, cbar_kwargs=None, **figure_kwargs):
    protein_states = data.columns.get_level_values(0).unique()

    n_subplots = len(protein_states)
    ncols = figure_kwargs.pop('ncols', min(cfg.getint('plotting', 'ncols'), n_subplots))
    nrows = figure_kwargs.pop('nrows', int(np.ceil(n_subplots / ncols)))
    figure_width = figure_kwargs.pop('width', cfg.getfloat('plotting', 'page_width')) / 25.4
    aspect = figure_kwargs.pop('aspect', cfg.getfloat('plotting', 'deltaG_aspect'))
    sharey = figure_kwargs.pop('sharey', 1)

    cmap_default, norm_default = get_cmap_norm_preset('vibrant', 10e3, 40e3)
    cmap = cmap or cmap_default
    cmap = pplt.Colormap(cmap)
    norm = norm or norm_default

    fig, axes = pplt.subplots(ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, sharey=sharey, **figure_kwargs)
    axes_iter = iter(axes)
    scatter_kwargs = scatter_kwargs or {}
    for state in protein_states:
        sub_df = data[state]
        ax = next(axes_iter)
        colorbar_scatter(ax, sub_df, cmap=cmap, norm=norm, cbar=False, **scatter_kwargs)

    for ax in axes_iter:
        ax.set_axis_off()

    # Set global ylims
    ylims = [lim for ax in axes if ax.axison for lim in ax.get_ylim()]
    axes.format(ylim=(np.max(ylims), np.min(ylims)), yticklabelloc='none', ytickloc='none')

    cbar_kwargs = cbar_kwargs or {}
    cbars = []
    cbar_norm = pplt.Norm('linear', norm.vmin*1e-3, norm.vmax*1e-3)
    for ax in axes:
        if not ax.axison:
            continue

        cbar = add_cbar(ax, cmap, cbar_norm, **cbar_kwargs)
        cbars.append(cbar)

    return fig, axes, cbars

#alias
deltaG_scatter_figure = dG_scatter_figure


def ddG_scatter_figure(data, reference=None, norm=None, cmap=None, scatter_kwargs=None, cbar_kwargs=None,
                       **figure_kwargs):
    protein_states = data.columns.get_level_values(0).unique()
    if reference is None:
        reference_state = protein_states[0]
    elif isinstance(reference, int):
        reference_state = protein_states[reference]
    elif reference in protein_states:
        reference_state = reference
    else:
        raise ValueError(f"Invalid value {reference!r} for 'reference'")

    dG_test = data.xs('deltaG', axis=1, level=1).drop(reference_state, axis=1)
    dG_ref = data[reference_state, 'deltaG']
    ddG = dG_test.subtract(dG_ref, axis=0)
    ddG.columns = pd.MultiIndex.from_product([ddG.columns, ['deltadeltaG']], names=['State', 'quantity'])

    cov_test = data.xs('covariance', axis=1, level=1).drop(reference_state, axis=1)**2
    cov_ref = data[reference_state, 'covariance']**2
    cov = cov_test.add(cov_ref, axis=0).pow(0.5)
    cov.columns = pd.MultiIndex.from_product([cov.columns, ['covariance']], names=['State', 'quantity'])

    combined = pd.concat([ddG, cov], axis=1)

    n_subplots = len(protein_states) - 1
    ncols = figure_kwargs.pop('ncols', min(cfg.getint('plotting', 'ncols'), n_subplots))
    nrows = figure_kwargs.pop('nrows', int(np.ceil(n_subplots / ncols)))
    figure_width = figure_kwargs.pop('width', cfg.getfloat('plotting', 'page_width')) / 25.4
    aspect = figure_kwargs.pop('aspect', cfg.getfloat('plotting', 'deltaG_aspect'))
    sharey = figure_kwargs.pop('sharey', 1)

    cmap_default, norm_default = get_cmap_norm_preset('PRGn', -10e3, 10e3)
    cmap = cmap or cmap_default
    cmap = pplt.Colormap(cmap)
    norm = norm or norm_default

    fig, axes = pplt.subplots(ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, sharey=sharey, **figure_kwargs)
    axes_iter = iter(axes)
    scatter_kwargs = scatter_kwargs or {}
    for state in protein_states:
        if state == reference_state:
            continue
        sub_df = combined[state]
        ax = next(axes_iter)
        colorbar_scatter(ax, sub_df, y='deltadeltaG', cmap=cmap, norm=norm, cbar=False, **scatter_kwargs)
        title = f'{state} - {reference_state}'
        ax.format(title=title)

    for ax in axes_iter:
        ax.set_axis_off()

    # Set global ylims
    ylim = np.abs([lim for ax in axes if ax.axison for lim in ax.get_ylim()]).max()
    axes.format(ylim=(ylim, -ylim), yticklabelloc='none', ytickloc='none')

    cbar_kwargs = cbar_kwargs or {}
    cbars = []
    cbar_norm = pplt.Norm('linear', norm.vmin*1e-3, norm.vmax*1e-3)
    for ax in axes:
        if not ax.axison:
            continue

        cbar = add_cbar(ax, cmap, cbar_norm, **cbar_kwargs)
        cbars.append(cbar)

    return fig, axes, cbars


deltadeltaG_scatter_figure = ddG_scatter_figure


def linear_bars(data, reference=None, field='deltaG', norm=None, cmap=None, labels=None, **figure_kwargs):
    #todo add sorting
    protein_states = data.columns.get_level_values(0).unique()

    if isinstance(reference, int):
        reference_state = protein_states[reference]
    elif reference in protein_states:
        reference_state = reference
    elif reference is None:
        reference_state = None
    else:
        raise ValueError(f"Invalid value {reference!r} for 'reference'")

    if reference_state:
        test = data.xs(field, axis=1, level=1).drop(reference_state, axis=1)
        ref = data[reference_state, field]
        plot_data = test.subtract(ref, axis=0)
        plot_data.columns = pd.MultiIndex.from_product([plot_data.columns, [field]], names=['State', 'quantity'])

        cmap_default, norm_default = get_cmap_norm_preset('PRGn', -10e3, 10e3)
        n_subplots = len(protein_states) - 1
    else:
        plot_data = data
        cmap_default, norm_default = get_cmap_norm_preset('vibrant', 10e3, 40e3)
        n_subplots = len(protein_states)

    cmap = cmap or cmap_default
    norm = norm or norm_default

    ncols = 1
    nrows = n_subplots
    figure_width = figure_kwargs.pop('width', cfg.getfloat('plotting', 'page_width')) / 25.4
    aspect = figure_kwargs.pop('aspect', cfg.getfloat('plotting', 'linear_bars_aspect'))
    cbar_width = figure_kwargs.pop('cbar_width', cfg.getfloat('plotting', 'cbar_width')) / 25.4

    fig, axes = pplt.subplots(nrows=nrows, ncols=ncols, aspect=aspect, width=figure_width, hspace=0)
    axes_iter = iter(axes)
    labels = labels or protein_states
    if len(labels) != len(protein_states):
        raise ValueError('Number of labels provided must be equal to the number of protein states')
    for label, state in zip(labels, protein_states):
        if state == reference_state:
            continue

        values = plot_data[state, field]
        rmin, rmax = values.index.min(), values.index.max()
        extent = [rmin - 0.5, rmax + 0.5, 0, 1]

        img = np.expand_dims(values, 0)

        ax = next(axes_iter)
        from matplotlib.axes import Axes
        Axes.imshow(ax, norm(img), aspect='auto', cmap=cmap, vmin=0, vmax=1, interpolation='None',
                    extent=extent)

        # ax.imshow(img, aspect='auto', cmap=cmap, norm=norm, interpolation='None', discrete=False,
        #             extent=extent)
        ax.format(yticks=[])
        ax.text(1.02, 0.5, label, horizontalalignment='left',
                  verticalalignment='center', transform=ax.transAxes)

    axes.format(xlabel=r_xlabel)

    sclf = 1e-3 # todo kwargs / check value of filed
    cmap_norm = copy(norm)
    cmap_norm.vmin *= sclf
    cmap_norm.vmax *= sclf

    if field == 'deltaG':
        label = dG_ylabel
    elif field == 'deltaG' and reference_state:
        label = ddG_ylabel
    else:
        label = ''

    fig.colorbar(cmap, norm=cmap_norm, loc='b', label=label, width=cbar_width)

    return fig, axes


def rainbowclouds(data, reference=None, field='deltaG', norm=None, cmap=None, update_rc=True, **figure_kwargs):
    # todo add sorting
    if update_rc:
        plt.rcParams["image.composite_image"] = False

    protein_states = data.columns.get_level_values(0).unique()

    if isinstance(reference, int):
        reference_state = protein_states[reference]
    elif reference in protein_states:
        reference_state = reference
    elif reference is None:
        reference_state = None
    else:
        raise ValueError(f"Invalid value {reference!r} for 'reference'")

    if reference_state:
        test = data.xs(field, axis=1, level=1).drop(reference_state, axis=1)
        ref = data[reference_state, field]
        plot_data = test.subtract(ref, axis=0)
        plot_data.columns = pd.MultiIndex.from_product([plot_data.columns, [field]], names=['State', 'quantity'])

        cmap_default, norm_default = get_cmap_norm_preset('PRGn', -10e3, 10e3)
    else:
        plot_data = data
        cmap_default, norm_default = get_cmap_norm_preset('vibrant', 10e3, 40e3)


    cmap = cmap or cmap_default
    norm = norm or norm_default
    plot_data = plot_data.xs(field, axis=1, level=1)

    #scaling
    plot_data *= 1e-3
    norm.vmin = norm.vmin * 1e-3
    norm.vmax = norm.vmax * 1e-3

    f_data = [plot_data[column].dropna().to_numpy() for column in plot_data.columns]  # todo make funcs accept dataframes
    f_labels = plot_data.columns

    ncols = 1
    nrows = 1
    figure_width = figure_kwargs.pop('width', cfg.getfloat('plotting', 'page_width')) / 25.4
    aspect = figure_kwargs.pop('aspect', cfg.getfloat('plotting', 'rainbow_aspect'))

    boxplot_width = 0.1
    orientation = 'vertical'

    strip_kwargs = dict(offset=0.0, orientation=orientation, s=2, colors='k', jitter=0.2, alpha=0.25)
    kde_kwargs = dict(linecolor='k', offset=0.15, orientation=orientation, fillcolor=False, fill_cmap=cmap,
                      fill_norm=norm, y_scale=None, y_norm=0.4, linewidth=1)
    boxplot_kwargs = dict(offset=0.2, sym='', linewidth=1., linecolor='k', orientation=orientation,
                          widths=boxplot_width)

    fig, axes = pplt.subplots(nrows=nrows, ncols=ncols, width=figure_width, aspect=aspect, hspace=0)
    ax = axes[0]
    stripplot(f_data, ax=ax, **strip_kwargs)
    kdeplot(f_data, ax=ax, **kde_kwargs)
    boxplot(f_data, ax=ax, **boxplot_kwargs)
    label_axes(f_labels, ax=ax, rotation=45)
    if field == 'deltaG':
        label = dG_ylabel
    elif field == 'deltaG' and reference_state:
        label = ddG_ylabel
    else:
        label = ''
    ax.format(xlim=(-0.75, len(f_data) - 0.5), ylabel=label, yticklabelloc='left', ytickloc='left',
              ylim=ax.get_ylim()[::-1])

    add_cbar(ax, cmap, norm)

    return fig, ax


def colorbar_scatter(ax, data, y='deltaG', yerr='covariance', cmap=None, norm=None, cbar=True, **kwargs):
    #todo refactor to colorbar_scatter?
    #todo custom ylims? scaling?
    if y == 'deltaG':
        cmap_default, norm_default = get_cmap_norm_preset('vibrant', 10e3, 40e3)
        sclf = 1e-3  # deltaG are given in J/mol but plotted in kJ/mol
    elif y == 'deltadeltaG':
        cmap_default, norm_default = get_cmap_norm_preset('PRGn', -10e3, 10e3)
        sclf = 1e-3
    else:
        if cmap is None or norm is None:
            raise ValueError("No valid `cmap` or `norm` is given.")
        sclf = 1e-3


    cmap = cmap or cmap_default
    cmap = pplt.Colormap(cmap)
    norm = norm or norm_default

    colors = cmap(norm(data[y]))

    #todo errorbars using proplot kwargs?
    errorbar_kwargs = {**ERRORBAR_KWARGS, **kwargs.pop('errorbar_kwargs', {})}
    scatter_kwargs = {**SCATTER_KWARGS, **kwargs}
    ax.scatter(data.index, data[y]*sclf, color=colors, **scatter_kwargs)
    with autoscale_turned_off(ax):
        ax.errorbar(data.index, data[y]*sclf, yerr=data[yerr]*sclf, zorder=-1,
                    **errorbar_kwargs)
    ax.set_xlabel(r_xlabel)
    # Default y labels
    labels = {'deltaG': dG_ylabel, 'deltadeltaG': ddG_ylabel}
    label = labels.get(y, '')
    ax.set_ylabel(label)

    ylim = ax.get_ylim()
    if (ylim[0] < ylim[1]) and y == 'deltaG':
        ax.set_ylim(*ylim[::-1])
    elif y == 'deltadeltaG':
        ylim = np.max(np.abs(ylim))
        ax.set_ylim(ylim, -ylim)


    if cbar:
        cbar_norm = copy(norm)
        cbar_norm.vmin *= sclf
        cbar_norm.vmax *= sclf
        cbar = add_cbar(ax, cmap, cbar_norm)
    else:
        cbar = None

    return cbar


def cmap_norm_from_nodes(colors, nodes, bad=None):
    nodes = np.array(nodes)
    if not np.all(np.diff(nodes) > 0):
        raise ValueError("Node values must be monotonically increasing")

    norm = pplt.Norm('linear', vmin=nodes.min(), vmax=nodes.max(), clip=True)
    color_spec = list(zip(norm(nodes), colors))
    cmap = pplt.Colormap(color_spec)
    bad = bad or cfg.get('plotting', 'no_coverage')
    cmap.set_bad(bad)

    return cmap, norm


def get_cmap_norm_preset(name, vmin, vmax):
    # Paul Tol colour schemes: https://personal.sron.nl/~pault/#sec:qualitative

    #todo warn if users use diverging colors with non diverging vmin/vmax?
    colors, bad = get_color_scheme(name)
    nodes = np.linspace(vmin, vmax, num=len(colors), endpoint=True)

    cmap, norm = cmap_norm_from_nodes(colors, nodes, bad)

    return cmap, norm


def get_color_scheme(name):
    # Paul Tol colour schemes: https://personal.sron.nl/~pault/#sec:qualitative
    if name == 'rgb':
        colors = ['#0000ff', '#00ff00', '#ff0000']  # red, green, blue
        bad = '#8c8c8c'
    elif name == 'bright':
        colors = ['#ee6677', '#288833', '#4477aa']
        bad = '#bbbbbb'
    elif name == 'vibrant':
        colors = ['#CC3311', '#009988', '#0077BB']
        bad = '#bbbbbb'
    elif name == 'muted':
        colors = ['#882255', '#117733', '#332288']
        bad = '#dddddd'
    elif name == 'pale':
        colors = ['#ffcccc', '#ccddaa', '#bbccee']
        bad = '#dddddd'
    elif name == 'dark':
        colors = ['#663333', '#225522', '#222255']
        bad = '#555555'
    elif name == 'delta':  # Original ddG colors
        colors = ['#006d2c', '#ffffff', '#54278e']  # Green, white, purple (flexible, no change, rigid)
        bad = '#ffee99'
    elif name == 'sunset':
        colors = ['#a50026', '#dd3d2d', '#f67e4b', '#fdb366', '#feda8b', '#eaeccc', '#c2e4ef', '#98cae1', '#6ea6cd',
                  '#4a7bb7', '#364b9a']
        bad = '#ffffff'
    elif name == 'BuRd':
        colors = ['#b2182b', '#d6604d', '#f4a582', '#fddbc7', '#f7f7f7', '#d1e5f0', '#92c5de', '#4393c3', '#2166ac']
        bad = '#ffee99'
    elif name == 'PRGn':
        colors = ['#1b7837', '#5aae61', '#acd39e', '#d9f0d3', '#f7f7f7', '#e7d4e8', '#c2a5cf', '#9970ab', '#762a83']
        bad = '#ffee99'
    else:
        raise ValueError(f"Color scheme '{name}' not found")

    return colors, bad


def pymol_figures(data, output_path, pdb_file, reference=None, field='deltaG', cmap=None, norm=None, extent=None,
                  orient=True, views=None, name_suffix='',
                  additional_views=None, img_size=(640, 640)):

    protein_states = data.columns.get_level_values(0).unique()

    if isinstance(reference, int):
        reference_state = protein_states[reference]
    elif reference in protein_states:
        reference_state = reference
    elif reference is None:
        reference_state = None
    else:
        raise ValueError(f"Invalid value {reference!r} for 'reference'")

    if reference_state:
        test = data.xs(field, axis=1, level=1).drop(reference_state, axis=1)
        ref = data[reference_state, field]
        plot_data = test.subtract(ref, axis=0)
        plot_data.columns = pd.MultiIndex.from_product([plot_data.columns, [field]], names=['State', 'quantity'])

        cmap_default, norm_default = get_cmap_norm_preset('PRGn', -10e3, 10e3)
    else:
        plot_data = data
        cmap_default, norm_default = get_cmap_norm_preset('vibrant', 10e3, 40e3)

    cmap = cmap or cmap_default
    norm = norm or norm_default
    #plot_data = plot_data.xs(field, axis=1, level=1)

    for state in protein_states:
        if state == reference_state:
            continue

        values = plot_data[state, field]
        rmin, rmax = extent or [None, None]
        rmin = rmin or values.index.min()
        rmax = rmax or values.index.max()

        values = values.reindex(pd.RangeIndex(rmin, rmax+1, name='r_number'))
        colors = apply_cmap(values, cmap, norm)
        name = f'pymol_ddG_{state}' if reference_state else f'pymol_dG_{state}'
        name += name_suffix
        pymol_render(output_path, pdb_file, colors, name=name, orient=orient, views=views, additional_views=additional_views,
                     img_size=img_size)


def pymol_render(output_path, pdb_file, colors, name='Pymol render', orient=True, views=None, additional_views=None, img_size=(640, 640)):
    if cmd is None:
        raise ModuleNotFoundError("Pymol module is not installed")

    px, py = img_size

    cmd.reinitialize()
    cmd.load(pdb_file)
    if orient:
        cmd.orient()
    cmd.set('antialias', 2)
    cmd.set('fog', 0)

    color_pymol(colors, cmd)

    if views:
        for i, view in enumerate(views):
            cmd.set_view(view)
            cmd.ray(px, py, renderer=0, antialias=2)
            output_file = output_path / f'{name}_view_{i}.png'
            cmd.png(str(output_file))

    else:
        cmd.ray(px, py, renderer=0, antialias=2)
        output_file = output_path / f'{name}_xy.png'
        cmd.png(str(output_file))

        cmd.rotate('x', 90)

        cmd.ray(px, py, renderer=0, antialias=2)
        output_file = output_path / f'{name}_xz.png'
        cmd.png(str(output_file))

        cmd.rotate('z', -90)

        cmd.ray(px, py, renderer=0, antialias=2)
        output_file = output_path / f'{name}_yz.png'
        cmd.png(str(output_file))

        additional_views = additional_views or []

        for i, view in enumerate(additional_views):
            cmd.set_view(view)
            cmd.ray(px, py, renderer=0, antialias=2)
            output_file = output_path / f'{name}_view_{i}.png'
            cmd.png(str(output_file))


def add_cbar(ax, cmap, norm, **kwargs):
    """Truncate or append cmap such that it covers axes limit and and colorbar to axes"""

    cmap = pplt.Colormap(cmap)

    vmin, vmax = norm.vmin, norm.vmax
    ylim = ax.get_ylim()
    ymin, ymax = np.min(ylim), np.max(ylim)

    nodes = [ymin, vmin, vmax, ymax]
    all_ratios = np.diff(nodes)
    idx = np.nonzero(all_ratios > 0)
    all_cmaps = np.array([pplt.Colormap([cmap(0.)]), cmap, pplt.Colormap([cmap(1.)])])
    cmaps = all_cmaps[idx]
    ratios = all_ratios[idx]
    if len(cmaps) >= 2:
        new_cmap = cmaps[0].append(*cmaps[1:], ratios=ratios)
    else:
        new_cmap = cmap
    reverse = ylim[0] > ylim[1]

    new_total_length = np.sum(ratios)
    left = np.max([-all_ratios[0] / new_total_length, 0.])
    right = np.min([1 + all_ratios[-1] / new_total_length, 1.])

    new_cmap = new_cmap.truncate(left=left, right=right)
    new_norm = pplt.Norm('linear', vmin=ymin, vmax=ymax)

    cbar_kwargs = {**CBAR_KWARGS, **kwargs}
    cbar = ax.colorbar(new_cmap, norm=new_norm, reverse=reverse, **cbar_kwargs)

    return cbar


#https://stackoverflow.com/questions/38629830/how-to-turn-off-autoscaling-in-matplotlib-pyplot
@contextmanager
def autoscale_turned_off(ax=None):
  ax = ax or plt.gca()
  lims = [ax.get_xlim(), ax.get_ylim()]
  yield
  ax.set_xlim(*lims[0])
  ax.set_ylim(*lims[1])


def stripplot(data, ax=None, jitter=0.25, colors=None, offset=0., orientation='vertical', **scatter_kwargs):
    ax = ax or plt.gca()
    color_list = _prepare_colors(colors, len(data))

    for i, (d, color) in enumerate(zip(data, color_list)):
        jitter_offsets = (np.random.rand(d.size) - 0.5) * jitter
        cat_var = i * np.ones_like(d) + jitter_offsets + offset  # categorical axis variable
        if orientation == 'vertical':
            ax.scatter(cat_var, d, color=color, **scatter_kwargs)
        elif orientation == 'horizontal':
            ax.scatter(d, len(data) - cat_var, color=color, **scatter_kwargs)


def _prepare_colors(colors, N):
    if not isinstance(colors, list):
        return [colors]*N
    else:
        return colors


# From joyplot
def _x_range(data, extra=0.2):
    """ Compute the x_range, i.e., the values for which the
        density will be computed. It should be slightly larger than
        the max and min so that the plot actually reaches 0, and
        also has a bit of a tail on both sides.
    """
    try:
        sample_range = np.nanmax(data) - np.nanmin(data)
    except ValueError:
        return []
    if sample_range < 1e-6:
        return [np.nanmin(data), np.nanmax(data)]
    return np.linspace(np.nanmin(data) - extra*sample_range,
                       np.nanmax(data) + extra*sample_range, 1000)


def kdeplot(data, ax=None, offset=0., orientation='vertical',
            linecolor=None, linewidth=None, zero_line=True, x_extend=1e-3, y_scale=None, y_norm=None, fillcolor=False, fill_cmap=None,
            fill_norm=None):
    assert not (y_scale and y_norm), "Cannot set both 'y_scale' and 'y_norm'"
    y_scale = 1. if y_scale is None else y_scale

    color_list = _prepare_colors(linecolor, len(data))

    for i, (d, color) in enumerate(zip(data, color_list)):
        #todo remove NaNs?

        # Perhaps also borrow this part from joyplot
        kde_func = kde.gaussian_kde(d)
        kde_x = _x_range(d, extra=0.4)
        kde_y = kde_func(kde_x)*y_scale
        if y_norm:
            kde_y = y_norm*kde_y / kde_y.max()
        bools = kde_y > x_extend * kde_y.max()
        kde_x = kde_x[bools]
        kde_y = kde_y[bools]

        cat_var = len(data) - i + kde_y + offset # x in horizontal
        cat_var_zero = (len(data) - i)*np.ones_like(kde_y) + offset

        # x = i * np.ones_like(d) + jitter_offsets + offset  # 'x' like, could be y axis
        if orientation == 'horizontal':
            plot_x = kde_x
            plot_y = cat_var
            img_data = kde_x.reshape(1, -1)
        elif orientation == 'vertical':
            plot_x = len(data) - cat_var
            plot_y = kde_x
            img_data = kde_x[::-1].reshape(-1, 1)
        else:
            raise ValueError(f"Invalid value '{orientation}' for 'orientation'")

        line, = ax.plot(plot_x, plot_y, color=color, linewidth=linewidth)
        if zero_line:
            ax.plot([plot_x[0], plot_x[-1]], [plot_y[0], plot_y[-1]], color=line.get_color(), linewidth=linewidth)

        if fillcolor:
            #todo refactor to one if/else orientation
            color = line.get_color() if fillcolor is True else fillcolor
            if orientation == 'horizontal':
                ax.fill_between(kde_x, plot_y, np.linspace(plot_y[0], plot_y[-1], num=plot_y.size, endpoint=True),
                                color=color)
            elif orientation == 'vertical':
                ax.fill_betweenx(kde_x, len(data) - cat_var, len(data) - cat_var_zero, color=color)

        if fill_cmap:
            fill_norm = fill_norm or (lambda x: x)
            color_img = fill_norm(img_data)

            xmin, xmax = np.min(plot_x), np.max(plot_x)
            ymin, ymax = np.min(plot_y), np.max(plot_y)
            extent = [xmin-offset, xmax-offset, ymin, ymax] if orientation == 'horizontal' else [xmin, xmax, ymin-offset, ymax-offset]
            im = Axes.imshow(ax, color_img, aspect='auto', cmap=fill_cmap, extent=extent)  # left, right, bottom, top
            fill_line, = ax.fill(plot_x, plot_y, facecolor='none')
            im.set_clip_path(fill_line)


def boxplot(data, ax, offset=0., orientation='vertical', widths=0.25, linewidth=None, linecolor=None, **kwargs):
    if orientation == 'vertical':
        vert = True
        positions = np.arange(len(data)) + offset
    elif orientation == 'horizontal':
        vert = False
        positions = len(data) - np.arange(len(data)) - offset
    else:
        raise ValueError(f"Invalid value '{orientation}' for 'orientation', options are 'horizontal' or 'vertical'")

    #todo for loop
    boxprops = kwargs.pop('boxprops', {})
    whiskerprops = kwargs.pop('whiskerprops', {})
    medianprops = kwargs.pop('whiskerprops', {})

    boxprops['linewidth'] = linewidth
    whiskerprops['linewidth'] = linewidth
    medianprops['linewidth'] = linewidth

    boxprops['color'] = linecolor
    whiskerprops['color'] = linecolor
    medianprops['color'] = linecolor

    Axes.boxplot(ax, data, vert=vert, positions=positions, widths=widths, boxprops=boxprops, whiskerprops=whiskerprops,
               medianprops=medianprops, **kwargs)


def label_axes(labels, ax, offset=0., orientation='vertical', **kwargs):
    #todo check offset sign
    if orientation == 'vertical':
        ax.set_xticks(np.arange(len(labels)) + offset)
        ax.set_xticklabels(labels, **kwargs)
    elif orientation == 'horizontal':
        ax.set_yticks(len(labels) - np.arange(len(labels)) + offset)
        ax.set_yticklabels(labels, **kwargs)


def plot_fitresults(fitresult_path, reference=None, plots='all', renew=False, cmap_and_norm=None, output_path=None,
                    output_type='.png', **save_kwargs):
    """

    Parameters
    ----------
    fitresult_path
    plots
    renew
    cmap_and_norm: :obj:`dict`, optional
        Dictionary with cmap and norms to use. If `None`, reverts to defaults.
        Dict format: {'dG': (cmap, norm), 'ddG': (cmap, norm)}

    output_type: list or str

    Returns
    -------

    """
    # batch results only
    history_path = fitresult_path / 'model_history.csv'
    output_path = output_path or fitresult_path
    output_type = list([output_type]) if isinstance(output_type, str) else output_type
    fitresult = load_fitresult(fitresult_path)

    protein_states = fitresult.output.df.columns.get_level_values(0).unique()

    if isinstance(reference, int):
        reference_state = protein_states[reference]
    elif reference in protein_states:
        reference_state = reference
    elif reference is None:
        reference_state = None
    else:
        raise ValueError(f"Invalid value {reference!r} for 'reference'")

    cmap_and_norm = cmap_and_norm or {}
    dG_cmap, dG_norm = cmap_and_norm.get('dG', (None, None))
    dG_cmap_default, dG_norm_default = get_cmap_norm_preset('vibrant', 10e3, 40e3)
    ddG_cmap, ddG_norm = cmap_and_norm.get('ddG', (None, None))
    ddG_cmap_default, ddG_norm_default = get_cmap_norm_preset('PRGn', -10e3, 10e3)
    dG_cmap = ddG_cmap or dG_cmap_default
    dG_norm = dG_norm or dG_norm_default
    ddG_cmap = ddG_cmap or ddG_cmap_default
    ddG_norm = ddG_norm or ddG_norm_default

    check_exists = lambda x: False if renew else x.exists()


    if plots == 'all':
        plots = ['loss', 'rfu_coverage', 'rfu_scatter', 'dG_scatter', 'ddG_scatter', 'linear_bars', 'rainbowclouds']


    # def check_update(pth, fname, extensions, renew):
    #     # Returns True if the target graph should be renewed or not
    #     if renew:
    #         return True
    #     else:
    #         pths = [pth / (fname + ext) for ext in extensions]
    #         return any([not pth.exists() for pth in pths])

    # plots = [p for p in plots if check_update(output_path, p, output_type, renew)]

    if 'loss' in plots:
        loss_df = fitresult.losses
        loss_df.plot()

        mse_loss = loss_df['mse_loss']
        reg_loss = loss_df.iloc[:, 1:].sum(axis=1)
        reg_percentage = 100*reg_loss / (mse_loss + reg_loss)
        fig = plt.gcf()
        ax = plt.gca()
        ax1 = ax.twinx()
        reg_percentage.plot(ax=ax1, color='k')
        ax1.set_xlim(0, None)
        for ext in output_type:
            f_out = output_path / ('loss' + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if 'rfu_coverage' in plots:
        for hdxm in fitresult.data_obj:
            fig, axes, cbar_ax = peptide_coverage_figure(hdxm.data)
            for ext in output_type:
                f_out = output_path / (f'rfu_coverage_{hdxm.name}' + ext)
                plt.savefig(f_out)
            plt.close(fig)

    #todo rfu_scatter_timepoint

    if 'rfu_scatter' in plots:
        fig, axes, cbar = residue_scatter_figure(fitresult.data_obj)
        for ext in output_type:
            f_out = output_path / (f'rfu_scatter' + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if 'dG_scatter' in plots:
        fig, axes, cbars = dG_scatter_figure(fitresult.output.df, cmap=dG_cmap, norm=dG_norm)
        for ext in output_type:
            f_out = output_path / (f'dG_scatter' + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if 'ddG_scatter' in plots:
        fig, axes, cbars = ddG_scatter_figure(fitresult.output.df, reference=reference, cmap=ddG_cmap, norm=ddG_norm)
        for ext in output_type:
            f_out = output_path / (f'ddG_scatter' + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if 'linear_bars' in plots:
        fig, axes = linear_bars(fitresult.output.df)
        for ext in output_type:
            f_out = output_path / (f'dG_linear_bars' + ext)
            plt.savefig(f_out)
        plt.close(fig)

        if reference_state:
            fig, axes = linear_bars(fitresult.output.df, reference=reference)
            for ext in output_type:
                f_out = output_path / (f'ddG_linear_bars' + ext)
                plt.savefig(f_out)
            plt.close(fig)

    if 'rainbowclouds' in plots:
        fig, ax = rainbowclouds(fitresult.output.df)
        for ext in output_type:
            f_out = output_path / (f'dG_rainbowclouds' + ext)
            plt.savefig(f_out)
        plt.close(fig)

        if reference_state:
            fig, axes = rainbowclouds(fitresult.output.df, reference=reference)
            for ext in output_type:
                f_out = output_path / (f'ddG_rainbowclouds' + ext)
                plt.savefig(f_out)
            plt.close(fig)


    #
    # if 'history' in plots:
    #     for h_df, name in zip(history_list, names):
    #         output_path = fitresult_path / f'{name}history.png'
    #         if check_exists(output_path):
    #             break
    #
    #         num = len(h_df.columns)
    #         max_epochs = max([int(c) for c in h_df.columns])
    #
    #         cmap = mpl.cm.get_cmap('winter')
    #         norm = mpl.colors.Normalize(vmin=1, vmax=max_epochs)
    #         colors = iter(cmap(np.linspace(0, 1, num=num)))
    #
    #         fig, axes = pplt.subplots(nrows=1, width=width, aspect=aspect)
    #         ax = axes[0]
    #         for key in h_df:
    #             c = next(colors)
    #             to_hex(c)
    #
    #             ax.scatter(h_df.index, h_df[key] * 1e-3, color=to_hex(c), **scatter_kwargs)
    #         ax.format(xlabel=r_xlabel, ylabel=dG_ylabel)
    #
    #         values = np.linspace(0, max_epochs, endpoint=True, num=num)
    #         colors = cmap(norm(values))
    #         tick_labels = np.linspace(0, max_epochs, num=5)
    #
    #         cbar = fig.colorbar(colors, values=values, ticks=tick_labels, space=0, width=cbar_width, label='Epochs')
    #         ax.format(yticklabelloc='None', ytickloc='None')
    #
    #         plt.savefig(output_path)
    #         plt.close(fig)


