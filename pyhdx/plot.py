from contextlib import contextmanager
from copy import copy
from pathlib import Path
from typing import Type, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import proplot as pplt
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from scipy.stats import kde
from tqdm import tqdm
import colorcet as cc

from pyhdx.config import cfg
from pyhdx.fileIO import load_fitresult
from pyhdx.support import (
    autowrap,
    color_pymol,
    apply_cmap,
    multiindex_astype,
    multiindex_set_categories,
)
from pyhdx.tol_colors import tol_cmap

try:
    from pymol import cmd
except ModuleNotFoundError:
    cmd = None

dG_ylabel = "ΔG (kJ/mol)"
ddG_ylabel = "ΔΔG (kJ/mol)"
r_xlabel = "Residue Number"
NO_COVERAGE = "#8c8c8c"

ERRORBAR_KWARGS = {
    "fmt": "o",
    "ecolor": "k",
    "elinewidth": 0.3,
    "markersize": 0,
    "alpha": 0.75,
    "capthick": 0.3,
    "capsize": 0.0,
}

SCATTER_KWARGS = {"s": 7}

RECT_KWARGS = {"linewidth": 0.5, "linestyle": "-", "edgecolor": "k"}

CBAR_KWARGS = {
    "space": 0,
    "width": cfg.getfloat("plotting", "cbar_width") / 25.4,
    "tickminor": True,
}


def peptide_coverage_figure(
    data: pd.DataFrame,
    wrap: int = None,
    cmap: Union[pplt.Colormap, mpl.colors.Colormap, str, tuple, dict] = "turbo",
    norm: Type[mpl.colors.Normalize] = None,
    color_field: str = "rfu",
    subplot_field: str = "exposure",
    rect_fields: tuple = ("start", "end"),
    rect_kwargs: dict = None,
    **figure_kwargs,
) -> tuple:

    subplot_values = data[subplot_field].unique()
    sub_dfs = {
        value: data.query(f"`{subplot_field}` == {value}") for value in subplot_values
    }

    n_subplots = len(subplot_values)

    ncols = figure_kwargs.pop("ncols", min(cfg.getint("plotting", "ncols"), n_subplots))
    nrows = figure_kwargs.pop("nrows", int(np.ceil(n_subplots / ncols)))
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    cbar_width = (
        figure_kwargs.pop("cbar_width", cfg.getfloat("plotting", "cbar_width")) / 25.4
    )
    aspect = figure_kwargs.pop(
        "aspect", cfg.getfloat("plotting", "peptide_coverage_aspect")
    )

    cmap = pplt.Colormap(cmap)
    norm = norm or pplt.Norm("linear", vmin=0, vmax=1)

    start_field, end_field = rect_fields
    if wrap is None:
        wrap = max(
            [
                autowrap(sub_df[start_field], sub_df[end_field])
                for sub_df in sub_dfs.values()
            ]
        )

    # TODO refaspect for proplot >= 0.9.5
    fig, axes = pplt.subplots(
        ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, **figure_kwargs
    )
    rect_kwargs = rect_kwargs or {}
    axes_iter = iter(axes)
    for value, sub_df in sub_dfs.items():
        ax = next(axes_iter)
        peptide_coverage(
            ax,
            sub_df,
            cmap=cmap,
            norm=norm,
            color_field=color_field,
            wrap=wrap,
            cbar=False,
            **rect_kwargs,
        )
        ax.format(title=f"{subplot_field}: {value}")

    for ax in axes_iter:
        ax.axis("off")

    start, end = data[start_field].min(), data[end_field].max()
    pad = 0.05 * (end - start)
    axes.format(xlim=(start - pad, end + pad), xlabel=r_xlabel)

    if not cmap.monochrome:
        cbar_ax = fig.colorbar(cmap, norm, width=cbar_width)
        cbar_ax.set_label(color_field, labelpad=-0)
    else:
        cbar_ax = None

    return fig, axes, cbar_ax


def peptide_coverage(
    ax,
    data,
    wrap=None,
    cmap=None,
    norm=None,
    color_field="rfu",
    rect_fields=("start", "end"),
    labels=False,
    cbar=True,
    cbar_kwargs=None,
    **kwargs,
):
    start_field, end_field = rect_fields
    data = data.sort_values(by=[start_field, end_field])

    wrap = wrap or autowrap(data[start_field], data[end_field])
    rect_kwargs = {**RECT_KWARGS, **kwargs}

    cmap_default, norm_default = CMAP_NORM_DEFAULTS.get(color_field, (None, None))

    cmap = pplt.Colormap(cmap) if cmap is not None else cmap_default
    norm = norm or norm_default
    i = -1
    for p_num, idx in enumerate(data.index):
        elem = data.loc[idx]
        if i < -wrap:
            i = -1

        if color_field is None:
            color = cmap(0.5)
        else:
            color = cmap(norm(elem[color_field]))

        width = elem[end_field] - elem[start_field]
        rect = Rectangle(
            (elem[start_field] - 0.5, i), width, 1, facecolor=color, **rect_kwargs
        )
        ax.add_patch(rect)
        if labels:
            rx, ry = rect.get_xy()
            cy = ry
            cx = rx
            ax.annotate(
                str(p_num), (cx, cy), color="k", fontsize=6, va="bottom", ha="right"
            )
        i -= 1

    ax.set_ylim(-wrap, 0)
    start, end = data[start_field].min(), data[end_field].max()
    pad = 0.05 * (end - start)
    ax.set_xlim(start - pad, end + pad)
    ax.set_yticks([])

    _cbar_kwargs = cbar_kwargs or {}
    cbar_kwargs = {
        **CBAR_KWARGS,
        **_cbar_kwargs,
    }  # TODO py39 dict union: d |= e or d | e

    if cbar and color_field:
        cbar_ax = ax.colorbar(cmap, norm=norm, **cbar_kwargs)
        cbar_ax.set_label(color_field, labelpad=-0)
    else:
        cbar_ax = None

    return cbar_ax


def residue_time_scatter_figure(
    hdxm,
    field="rfu",
    cmap="turbo",
    norm=None,
    scatter_kwargs=None,
    cbar_kwargs=None,
    **figure_kwargs,
):
    """per-residue per-exposure values for field  `field` by weighted averaging"""

    n_subplots = hdxm.Nt
    ncols = figure_kwargs.pop("ncols", min(cfg.getint("plotting", "ncols"), n_subplots))
    nrows = figure_kwargs.pop("nrows", int(np.ceil(n_subplots / ncols)))
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    aspect = figure_kwargs.pop(
        "aspect", cfg.getfloat("plotting", "residue_scatter_aspect")
    )
    cbar_width = (
        figure_kwargs.pop("cbar_width", cfg.getfloat("plotting", "cbar_width")) / 25.4
    )

    cmap = pplt.Colormap(cmap)  # todo allow None as cmap
    norm = norm or pplt.Norm("linear", vmin=0, vmax=1)

    fig, axes = pplt.subplots(
        ncols=ncols,
        nrows=nrows,
        width=figure_width,
        aspect=aspect,
        sharey=4,
        **figure_kwargs,
    )
    scatter_kwargs = scatter_kwargs or {}
    axes_iter = iter(axes)
    for hdx_tp in hdxm:
        ax = next(axes_iter)
        residue_time_scatter(
            ax, hdx_tp, field=field, cmap=cmap, norm=norm, cbar=False, **scatter_kwargs
        )  # todo cbar kwargs? (check with other methods)
        ax.format(title=f"exposure: {hdx_tp.exposure:.1f}")

    for ax in axes_iter:
        ax.axis("off")

    axes.format(xlabel=r_xlabel, ylabel=field)

    cbar_kwargs = cbar_kwargs or {}
    cbars = []
    for ax in axes:
        if not ax.axison:
            continue

        cbar = add_cbar(ax, cmap, norm, **cbar_kwargs)
        cbars.append(cbar)

    return fig, axes, cbars


def residue_time_scatter(
    ax, hdx_tp, field="rfu", cmap="turbo", norm=None, cbar=True, **kwargs
):
    # update cmap, norm defaults
    cmap = pplt.Colormap(cmap)  # todo allow None as cmap
    norm = norm or pplt.Norm("linear", vmin=0, vmax=1)
    cbar_width = kwargs.pop("cbar_width", cfg.getfloat("plotting", "cbar_width")) / 25.4

    scatter_kwargs = {**SCATTER_KWARGS, **kwargs}
    values = hdx_tp.weighted_average(field)
    colors = cmap(norm(values))
    ax.scatter(values.index, values, c=colors, **scatter_kwargs)

    if not cmap.monochrome and cbar:
        add_cbar(ax, cmap, norm, width=cbar_width)


def residue_scatter_figure(
    hdxm_set,
    field="rfu",
    cmap="viridis",
    norm=None,
    scatter_kwargs=None,
    **figure_kwargs,
):
    n_subplots = hdxm_set.Ns
    ncols = figure_kwargs.pop("ncols", min(cfg.getint("plotting", "ncols"), n_subplots))
    nrows = figure_kwargs.pop(
        "nrows", int(np.ceil(n_subplots / ncols))
    )  # todo disallow setting rows
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    cbar_width = (
        figure_kwargs.pop("cbar_width", cfg.getfloat("plotting", "cbar_width")) / 25.4
    )
    aspect = figure_kwargs.pop(
        "aspect", cfg.getfloat("plotting", "residue_scatter_aspect")
    )

    cmap = pplt.Colormap(cmap)
    if norm is None:
        tps = np.unique(np.concatenate([hdxm.timepoints for hdxm in hdxm_set]))
        tps = tps[np.nonzero(tps)]
        norm = pplt.Norm("log", vmin=tps.min(), vmax=tps.max())
    else:
        tps = np.unique(np.concatenate([hdxm.timepoints for hdxm in hdxm_set]))

    fig, axes = pplt.subplots(
        ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, **figure_kwargs
    )
    axes_iter = iter(axes)
    scatter_kwargs = scatter_kwargs or {}
    for hdxm in hdxm_set:
        ax = next(axes_iter)
        residue_scatter(
            ax, hdxm, cmap=cmap, norm=norm, field=field, cbar=False, **scatter_kwargs
        )
        ax.format(title=f"{hdxm.name}")

    for ax in axes_iter:
        ax.axis("off")

    # todo function for this?
    locator = pplt.Locator(norm(tps))
    cbar_ax = fig.colorbar(cmap, width=cbar_width, ticks=locator)
    formatter = pplt.Formatter("simple", precision=1)
    cbar_ax.ax.set_yticklabels([formatter(t) for t in tps])
    cbar_ax.set_label("Exposure time (s)", labelpad=-0)

    axes.format(xlabel=r_xlabel)

    return fig, axes, cbar_ax


# todo allow colorbar_scatter to take rfus
def residue_scatter(
    ax, hdxm, field="rfu", cmap="viridis", norm=None, cbar=True, **kwargs
):
    cmap = pplt.Colormap(cmap)
    tps = hdxm.timepoints[np.nonzero(hdxm.timepoints)]
    norm = norm or pplt.Norm("log", tps.min(), tps.max())

    cbar_width = kwargs.pop("cbar_width", cfg.getfloat("plotting", "cbar_width")) / 25.4
    scatter_kwargs = {**SCATTER_KWARGS, **kwargs}
    for hdx_tp in hdxm:
        if isinstance(norm, mpl.colors.LogNorm) and hdx_tp.exposure == 0.0:
            continue
        values = hdx_tp.weighted_average(field)
        color = cmap(norm(hdx_tp.exposure))
        scatter_kwargs["color"] = color
        ax.scatter(values.index, values, **scatter_kwargs)

    if cbar:
        locator = pplt.Locator(norm(tps))
        cbar_ax = ax.colorbar(cmap, width=cbar_width, ticks=locator)
        formatter = pplt.Formatter("simple", precision=1)
        cbar_ax.ax.set_yticklabels([formatter(t) for t in tps])
        cbar_ax.set_label("Exposure time (s)", labelpad=-0)


def dG_scatter_figure(
    data, cmap=None, norm=None, scatter_kwargs=None, cbar_kwargs=None, **figure_kwargs
):
    protein_states = data.columns.get_level_values(0).unique()

    n_subplots = len(protein_states)
    ncols = figure_kwargs.pop("ncols", min(cfg.getint("plotting", "ncols"), n_subplots))
    nrows = figure_kwargs.pop("nrows", int(np.ceil(n_subplots / ncols)))
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    aspect = figure_kwargs.pop("aspect", cfg.getfloat("plotting", "dG_aspect"))
    sharey = figure_kwargs.pop("sharey", 1)

    cmap_default, norm_default = CMAP_NORM_DEFAULTS["dG"]
    cmap = cmap or cmap_default
    cmap = pplt.Colormap(cmap)
    norm = norm or norm_default

    fig, axes = pplt.subplots(
        ncols=ncols,
        nrows=nrows,
        width=figure_width,
        aspect=aspect,
        sharey=sharey,
        **figure_kwargs,
    )
    axes_iter = iter(axes)
    scatter_kwargs = scatter_kwargs or {}
    for state in protein_states:
        sub_df = data[state]
        ax = next(axes_iter)
        colorbar_scatter(
            ax,
            sub_df,
            cmap=cmap,
            norm=norm,
            cbar=False,
            sclf=1e-3,
            invert_yaxis=True,
            **scatter_kwargs,
        )
        ax.format(title=f"{state}")

    for ax in axes_iter:
        ax.set_axis_off()

    # Set global ylims
    ylims = [lim for ax in axes if ax.axison for lim in ax.get_ylim()]
    axes.format(
        ylim=(np.max(ylims), np.min(ylims)), yticklabelloc="none", ytickloc="none"
    )

    cbar_kwargs = cbar_kwargs or {}
    cbars = []
    cbar_norm = pplt.Norm("linear", norm.vmin * 1e-3, norm.vmax * 1e-3)
    for ax in axes:
        if not ax.axison:
            continue

        cbar = add_cbar(ax, cmap, cbar_norm, **cbar_kwargs)
        cbars.append(cbar)

    return fig, axes, cbars


def ddG_scatter_figure(
    data,
    reference=None,
    cmap=None,
    norm=None,
    scatter_kwargs=None,
    cbar_kwargs=None,
    **figure_kwargs,
):
    protein_states = data.columns.get_level_values(0).unique()
    if reference is None:
        reference_state = protein_states[0]
    elif isinstance(reference, int):
        reference_state = protein_states[reference]
    elif reference in protein_states:
        reference_state = reference
    else:
        raise ValueError(f"Invalid value {reference!r} for 'reference'")

    dG_test = data.xs("dG", axis=1, level=1).drop(reference_state, axis=1)
    dG_ref = data[reference_state, "dG"]
    ddG = dG_test.subtract(dG_ref, axis=0)
    ddG.columns = pd.MultiIndex.from_product(
        [ddG.columns, ["ddG"]], names=["State", "quantity"]
    )

    cov_test = data.xs("covariance", axis=1, level=1).drop(reference_state, axis=1) ** 2
    cov_ref = data[reference_state, "covariance"] ** 2
    cov = cov_test.add(cov_ref, axis=0).pow(0.5)
    cov.columns = pd.MultiIndex.from_product(
        [cov.columns, ["covariance"]], names=["State", "quantity"]
    )

    combined = pd.concat([ddG, cov], axis=1)

    n_subplots = len(protein_states) - 1
    ncols = figure_kwargs.pop("ncols", min(cfg.getint("plotting", "ncols"), n_subplots))
    nrows = figure_kwargs.pop("nrows", int(np.ceil(n_subplots / ncols)))
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    aspect = figure_kwargs.pop("aspect", cfg.getfloat("plotting", "dG_aspect"))
    sharey = figure_kwargs.pop("sharey", 1)

    cmap_default, norm_default = CMAP_NORM_DEFAULTS["ddG"]
    cmap = cmap or cmap_default
    cmap = pplt.Colormap(cmap)
    norm = norm or norm_default

    fig, axes = pplt.subplots(
        ncols=ncols,
        nrows=nrows,
        width=figure_width,
        aspect=aspect,
        sharey=sharey,
        **figure_kwargs,
    )
    axes_iter = iter(axes)
    scatter_kwargs = scatter_kwargs or {}
    for state in protein_states:
        if state == reference_state:
            continue
        sub_df = combined[state]
        ax = next(axes_iter)
        colorbar_scatter(
            ax,
            sub_df,
            y="ddG",
            cmap=cmap,
            norm=norm,
            cbar=False,
            sclf=1e-3,
            invert_yaxis=True,
            symmetric=True,
            **scatter_kwargs,
        )
        title = f"{state} - {reference_state}"
        ax.format(title=title)

    for ax in axes_iter:
        ax.set_axis_off()

    # Set global ylims
    ylim = np.abs([lim for ax in axes if ax.axison for lim in ax.get_ylim()]).max()
    axes.format(ylim=(ylim, -ylim), yticklabelloc="none", ytickloc="none")

    cbar_kwargs = cbar_kwargs or {}
    cbars = []
    cbar_norm = pplt.Norm("linear", norm.vmin * 1e-3, norm.vmax * 1e-3)
    for ax in axes:
        if not ax.axison:
            continue

        cbar = add_cbar(ax, cmap, cbar_norm, **cbar_kwargs)
        cbars.append(cbar)

    return fig, axes, cbars


def peptide_mse_figure(
    peptide_mse, cmap=None, norm=None, rect_kwargs=None, **figure_kwargs
):
    n_subplots = len(peptide_mse.columns.unique(level=0))
    ncols = figure_kwargs.pop("ncols", min(cfg.getint("plotting", "ncols"), n_subplots))
    nrows = figure_kwargs.pop("nrows", int(np.ceil(n_subplots / ncols)))
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    aspect = figure_kwargs.pop("aspect", cfg.getfloat("plotting", "peptide_mse_aspect"))

    cmap = cmap or CMAP_NORM_DEFAULTS["mse"][0]

    fig, axes = pplt.subplots(
        ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, **figure_kwargs
    )
    axes_iter = iter(axes)
    cbars = []
    rect_kwargs = rect_kwargs or {}
    states = peptide_mse.columns.unique(level=0)
    for state in states:
        sub_df = peptide_mse[state].dropna(how="all").convert_dtypes()

        ax = next(axes_iter)
        vmax = sub_df["peptide_mse"].max()

        # todo allow global norm by kwargs
        norm = norm or pplt.Norm("linear", vmin=0, vmax=vmax)
        # color bar per subplot as norm differs
        # todo perhaps unify color scale? -> when global norm, global cbar
        cbar_ax = peptide_coverage(
            ax, sub_df, color_field="peptide_mse", norm=norm, cmap=cmap, **rect_kwargs
        )
        cbar_ax.set_label("MSE")
        cbars.append(cbar_ax)
        ax.format(xlabel=r_xlabel, title=f"{state}")

    return fig, axes, cbars


def loss_figure(fit_result, **figure_kwargs):
    ncols = 1
    nrows = 1
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    aspect = figure_kwargs.pop(
        "aspect", cfg.getfloat("plotting", "loss_aspect")
    )  # todo loss aspect also in config?

    fig, ax = pplt.subplots(
        ncols=ncols, nrows=nrows, width=figure_width, aspect=aspect, **figure_kwargs
    )
    fit_result.losses.plot(ax=ax)
    # ax.plot(fit_result.losses, legend='t')  # altnernative proplot plotting

    # ox = ax.alty()
    # reg_loss = fit_result.losses.drop('mse_loss', axis=1)
    # total = fit_result.losses.sum(axis=1)
    # perc = reg_loss.divide(total, axis=0) * 100
    # perc.plot(ax=ox)  #todo formatting (perc as --, matching colors, legend)
    #

    ax.format(xlabel="Number of epochs", ylabel="Loss")

    return fig, ax


def linear_bars_figure(
    data,
    reference=None,
    groupby=None,
    field="dG",
    norm=None,
    cmap=None,
    **figure_kwargs,
):
    if reference is None and field == "dG":
        cmap_default, norm_default = CMAP_NORM_DEFAULTS["dG"]
        ylabel = dG_ylabel
        sclf = 1e-3
    elif reference is not None and field == "dG":
        cmap_default, norm_default = CMAP_NORM_DEFAULTS["ddG"]
        ylabel = ddG_ylabel
        sclf = 1e-3
    elif reference is None and field == "rfu":
        cmap_default, norm_default = CMAP_NORM_DEFAULTS["rfu"]
        ylabel = "RFU"
        sclf = 1.0
    elif reference is not None and field == "rfu":
        cmap_default, norm_default = CMAP_NORM_DEFAULTS["drfu"]
        ylabel = "ΔRFU"
        sclf = 1.0
    else:
        cmap_default, norm_default = None, None
        ylabel = ""
        sclf = 1.0

    cmap = cmap or cmap_default
    norm = norm or norm_default

    if cmap is None:
        raise ValueError("No valid Colormap found")
    if norm is None:
        raise ValueError("No valid Norm found")

    fig, axes, cbar = linear_bars(
        data,
        reference=reference,
        groupby=groupby,
        field=field,
        norm=norm,
        cmap=cmap,
        sclf=sclf,
        **figure_kwargs,
    )

    cbar.set_label(ylabel)

    return fig, axes, cbar


def linear_bars(
    data,
    reference=None,
    groupby=None,
    field="dG",
    norm=None,
    cmap=None,
    sclf=1.0,
    sort=False,
    **figure_kwargs,
):

    if data.columns.nlevels == 2:
        data = data.copy()
        columns = pd.MultiIndex.from_tuples(
            [("", *tup) for tup in data.columns], names=["group"] + data.columns.names
        )
        data.columns = columns

    data = data.xs(level=-1, key=field, drop_level=False, axis=1)

    groupby = groupby or data.columns.names[0]
    grp_level = data.columns.names.index(groupby)
    bars_level = 1 - grp_level

    if isinstance(reference, int):
        reference = data.columns.get_level_values(level=bars_level)[reference]
    if reference:
        ref_data = data.xs(key=reference, axis=1, level=bars_level)
        sub = data.subtract(ref_data, axis=1).reorder_levels(data.columns.names, axis=1)

        # fix column dtypes to preserve , reorder levels back to original order
        columns = multiindex_astype(sub.columns, grp_level, "category")
        categories = list(data.columns.unique(level=grp_level))
        columns = multiindex_set_categories(
            columns, grp_level, categories, ordered=True
        )
        sub.columns = columns
        sub = sub.sort_index(axis=1, level=grp_level)
        sub = sub.drop(axis=1, level=bars_level, labels=reference)

        plot_data = sub
    else:
        plot_data = data

    groups = plot_data.groupby(level=groupby, axis=1).groups
    hspace = [elem for v in groups.values() for elem in [0] * (len(v) - 1) + [None]][
        :-1
    ]

    ncols = 1
    nrows = len(hspace) + 1
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    aspect = figure_kwargs.pop("aspect", cfg.getfloat("plotting", "linear_bars_aspect"))
    cbar_width = (
        figure_kwargs.pop("cbar_width", cfg.getfloat("plotting", "cbar_width")) / 25.4
    )

    fig, axes = pplt.subplots(
        nrows=nrows, ncols=ncols, aspect=aspect, width=figure_width, hspace=hspace
    )
    axes_iter = iter(axes)

    if sort in ["ascending", "descending"]:
        srt_groups = plot_data.groupby(level=bars_level, axis=1).groups
        values = [plot_data.loc[:, grp].mean().mean() for grp in srt_groups.values()]
        idx = np.argsort(values)
        if sort == "descending":
            idx = idx[::-1]
    else:
        idx = None

    for grp_name, items in groups.items():
        items = items.values[idx] if idx is not None else items.values
        for i, item in enumerate(items):  # these items need to be sorted
            values = plot_data.xs(
                key=(*item[: plot_data.columns.nlevels - 1], field), axis=1
            )
            rmin, rmax = values.index.min(), values.index.max()
            extent = [rmin - 0.5, rmax + 0.5, 0, 1]
            img = np.expand_dims(values, 0)

            ax = next(axes_iter)
            label = item[bars_level]
            from matplotlib.axes import Axes

            Axes.imshow(  # TODO use proplot pcolor or related function
                ax,
                norm(img),
                aspect="auto",
                cmap=cmap,
                vmin=0,
                vmax=1,
                interpolation="None",
                extent=extent,
            )
            ax.format(yticks=[])
            ax.text(
                1.02,
                0.5,
                label,
                horizontalalignment="left",
                verticalalignment="center",
                transform=ax.transAxes,
            )
            if i == 0:
                ax.format(title=grp_name)

    axes.format(xlabel=r_xlabel)

    cmap_norm = copy(norm)
    cmap_norm.vmin *= sclf
    cmap_norm.vmax *= sclf

    cbar = fig.colorbar(cmap, norm=cmap_norm, loc="b", width=cbar_width)

    return fig, axes, cbar


def rainbowclouds_figure(
    data,
    reference=None,
    field="dG",
    norm=None,
    cmap=None,
    update_rc=True,
    **figure_kwargs,
):
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
        plot_data.columns = pd.MultiIndex.from_product(
            [plot_data.columns, [field]], names=["State", "quantity"]
        )

        cmap_default, norm_default = CMAP_NORM_DEFAULTS["ddG"]
        ylabel = ddG_ylabel
    else:
        plot_data = data
        cmap_default, norm_default = CMAP_NORM_DEFAULTS["dG"]
        ylabel = dG_ylabel

    cmap = cmap or cmap_default
    norm = norm or norm_default
    plot_data = plot_data.xs(field, axis=1, level=1)

    # scaling
    plot_data *= 1e-3
    norm = copy(norm)
    norm.vmin = norm.vmin * 1e-3
    norm.vmax = norm.vmax * 1e-3

    f_data = [
        plot_data[column].dropna().to_numpy() for column in plot_data.columns
    ]  # todo make funcs accept dataframes
    f_labels = plot_data.columns

    ncols = 1
    nrows = 1
    figure_width = (
        figure_kwargs.pop("width", cfg.getfloat("plotting", "page_width")) / 25.4
    )
    aspect = figure_kwargs.pop(
        "aspect", cfg.getfloat("plotting", "rainbowclouds_aspect")
    )

    fig, axes = pplt.subplots(
        nrows=nrows, ncols=ncols, width=figure_width, aspect=aspect, hspace=0
    )
    ax = axes[0]

    cbar = rainbowclouds(
        ax,
        f_data,
        f_labels,
        cmap=cmap,
        norm=norm,
        invert_yaxis=True,
        cbar_kwargs=CBAR_KWARGS,
    )
    cbar.set_label(ylabel)
    ax.format(ytickloc="none")

    return fig, ax, cbar


def rainbowclouds(
    ax,
    f_data,
    f_labels,
    cmap=None,
    norm=None,
    invert_yaxis=False,
    format_kwargs=None,
    cbar_kwargs=None,
    strip_kwargs=None,
    kde_kwargs=None,
    boxplot_kwargs=None,
):
    boxplot_width = 0.1
    orientation = "vertical"

    _strip_kwargs = dict(
        offset=0.0, orientation=orientation, s=2, colors="k", jitter=0.2, alpha=0.25
    )
    _kde_kwargs = dict(
        linecolor="k",
        offset=0.15,
        orientation=orientation,
        fillcolor=False,
        fill_cmap=cmap,
        fill_norm=norm,
        y_scale=None,
        y_norm=0.4,
        linewidth=1,
    )
    _boxplot_kwargs = dict(
        offset=0.2,
        sym="",
        linewidth=1.0,
        linecolor="k",
        orientation=orientation,
        widths=boxplot_width,
    )

    strip_kwargs = _strip_kwargs.update(strip_kwargs) if strip_kwargs else _strip_kwargs
    kde_kwargs = _kde_kwargs.update(strip_kwargs) if kde_kwargs else _kde_kwargs
    boxplot_kwargs = (
        _boxplot_kwargs.update(strip_kwargs) if boxplot_kwargs else _boxplot_kwargs
    )

    stripplot(f_data, ax=ax, **strip_kwargs)
    kdeplot(f_data, ax=ax, **kde_kwargs)
    boxplot(f_data, ax=ax, **boxplot_kwargs)
    label_axes(f_labels, ax=ax, rotation=45)

    ylim = ax.get_ylim()[::-1] if invert_yaxis else None
    _format_kwargs = dict(
        xlim=(-0.75, len(f_data) - 0.5),
        yticklabelloc="left",
        ytickloc="left",
        ylim=ylim,
    )
    format_kwargs = (
        _format_kwargs.update(format_kwargs) if format_kwargs else _format_kwargs
    )

    ax.format(**format_kwargs)

    if cmap is not None:
        cbar_kwargs = cbar_kwargs or {}
        cbar = add_cbar(ax, cmap, norm, **cbar_kwargs)
    else:
        cbar = None

    return cbar


# todo this should also take Pd.Series?
def colorbar_scatter(
    ax,
    data,
    y="dG",
    yerr="covariance",
    cmap=None,
    norm=None,
    cbar=True,
    cbar_kwargs=None,
    invert_yaxis=None,
    symmetric=False,
    sclf=1e-3,
    **kwargs,
):
    # todo make error bars optional
    # todo custom ylims? scaling?
    try:
        cmap_default, norm_default = CMAP_NORM_DEFAULTS[y]
    except KeyError:
        cmap_default, norm_default = None, None

    cmap = cmap or cmap_default
    cmap = pplt.Colormap(cmap) if isinstance(cmap, str) else cmap
    norm = norm or norm_default

    if cmap is None or norm is None:
        raise ValueError("No valid `cmap` or `norm` is given.")

    colors = cmap(norm(data[y]))

    # todo errorbars using proplot kwargs?
    errorbar_kwargs = {**ERRORBAR_KWARGS, **kwargs.pop("errorbar_kwargs", {})}
    scatter_kwargs = {**SCATTER_KWARGS, **kwargs}
    ax.scatter(data.index, data[y] * sclf, color=colors, **scatter_kwargs)
    with autoscale_turned_off(ax):
        ax.errorbar(
            data.index,
            data[y] * sclf,
            yerr=data[yerr] * sclf,
            zorder=-1,
            **errorbar_kwargs,
        )
    ax.set_xlabel(r_xlabel)
    # Default y labels
    labels = {"dG": dG_ylabel, "ddG": ddG_ylabel}
    label = labels.get(y, "")
    ax.set_ylabel(label)

    # todo this function should be more general and should take `invert_yaxis` and `symmetric` kwargs

    ylim = ax.get_ylim()
    if (ylim[0] < ylim[1]) and invert_yaxis and symmetric:
        ylim = np.max(np.abs(ylim))
        ax.set_ylim(ylim, -ylim)
    elif (ylim[0] < ylim[1]) and invert_yaxis:
        ax.set_ylim(*ylim[::-1])
    elif y == "ddG":
        ylim = np.max(np.abs(ylim))
        ax.set_ylim(ylim, -ylim)

    if cbar:
        cbar_norm = copy(norm)
        cbar_norm.vmin *= sclf
        cbar_norm.vmax *= sclf
        cbar_kwargs = cbar_kwargs or {}
        cbar = add_cbar(ax, cmap, cbar_norm, **cbar_kwargs)
    else:
        cbar = None

    return cbar


def redundancy(ax, hdxm, cmap=None, norm=None):
    cmap = cmap or CMAP_NORM_DEFAULTS.cmaps["redundancy"]
    norm = norm or CMAP_NORM_DEFAULTS.norms["redundancy"]

    redundancy = hdxm.coverage.X.sum(axis=0).astype(float)
    redundancy[redundancy == 0] = np.nan
    x = hdxm.coverage.r_number
    img = np.expand_dims(redundancy, 0)

    collection = ax.pcolormesh(
        pplt.edges(x), np.array([0, 1]), img, extend="max", cmap=cmap, norm=norm,
    )
    ax.format(yticks=[], title="Redundancy")
    return collection


def resolution(ax, hdxm, cmap=None, norm=None):
    cmap = cmap or CMAP_NORM_DEFAULTS.cmaps["resolution"]
    norm = norm or CMAP_NORM_DEFAULTS.norms["resolution"]

    resolution = np.repeat(hdxm.coverage.block_length, hdxm.coverage.block_length)
    resolution = resolution.astype(float)
    resolution[hdxm.coverage.X.sum(axis=0) == 0] = np.nan
    # TODO collection = single_linear_bar(ax, hdxm.coverage.r_number, resolution, cmap, norm)
    x = hdxm.coverage.r_number
    img = np.expand_dims(resolution, 0)

    collection = ax.pcolormesh(
        pplt.edges(x), np.array([0, 1]), img, extend="max", cmap=cmap, norm=norm,
    )
    ax.format(yticks=[], title="Resolution")
    return collection


def single_linear_bar(ax, x, z, cmap, norm, **kwargs):
    """makes a linear bar plot on supplied axis with values z and corresponding x values x"""

    if isinstance(z, pd.Series):
        z = z.to_numpy()
    elif isinstance(z, pd.DataFrame):
        assert len(z.columns) == 1, "Can only plot dataframes with 1 column"
        z = z.to_numpy().squeeze()

    img = np.expand_dims(z, 0)

    collection = ax.pcolormesh(
        pplt.edges(x), np.array([0, 1]), img, cmap=cmap, norm=norm, **kwargs
    )
    ax.format(yticks=[])

    return collection


def add_mse_panels(
    axes,
    fit_result,
    cmap=None,
    norm=None,
    panel_kwargs=None,
    fig=None,
    cbar=False,
    cbar_kwargs=None,
):
    """adds linear bar panels to axes showing some property (usually MSE)"""

    residue_mse = fit_result.get_residue_mse()
    vmax = residue_mse.to_numpy().max()

    cmap = cmap or CMAP_NORM_DEFAULTS.cmaps["mse"]
    norm = norm or pplt.Norm("linear", vmin=0, vmax=vmax)

    collections = []
    for hdxm, ax in zip(fit_result.hdxm_set, axes):
        panel_kwargs = panel_kwargs or {}
        loc = panel_kwargs.pop("loc", "b")
        panel_kwargs = {"width": "5mm", "space": 0, **panel_kwargs}

        pn = ax.panel_axes(loc, **panel_kwargs)
        residue_values = residue_mse[hdxm.name]

        collection = single_linear_bar(
            pn,
            hdxm.coverage.r_number,
            residue_values.to_numpy().squeeze(),
            cmap=cmap,
            norm=norm,
        )
        collections.append(collection)

    if cbar:
        if fig is None:
            raise ValueError(
                "Must pass 'fig' keyword argument to add a global colorbar"
            )
        cbar_kwargs = cbar_kwargs or {}
        cbar_kwargs = {
            "width": CBAR_KWARGS["width"],
            "loc": "b",
            "length": 0.3,
            **cbar_kwargs,
        }

        cbar = fig.colorbar(cmap, norm=norm, **cbar_kwargs)
    else:
        cbar = None

    return collections, cbar


def cmap_norm_from_nodes(colors, nodes, bad=None, under=None, over=None):
    nodes = np.array(nodes)
    if not np.all(np.diff(nodes) > 0):
        raise ValueError("Node values must be monotonically increasing")

    norm = pplt.Norm("linear", vmin=nodes.min(), vmax=nodes.max(), clip=True)
    color_spec = list(zip(norm(nodes), colors))
    cmap = pplt.Colormap(color_spec)

    if bad is not None:
        cmap.set_bad(bad)
    if under is not None:
        cmap.set_under(under)
    if over is not None:
        cmap.set_over(over)

    return cmap, norm


def set_bad(cmap, color=NO_COVERAGE):
    cmap.set_bad(color)
    return cmap


class ColorTransforms(object):
    def __init__(self):
        # Lily blue colormap
        foldedness_cmap, foldedness_cmap = cmap_norm_from_nodes(
            colors=["#ffffff", "#00ffff", "#008080", "#0075ea", "#000008"],
            nodes=[0.0, 0.25, 0.5, 0.75, 1.0],
            bad=NO_COVERAGE,
        )

        self.norms = {
            "dG": pplt.Norm("linear", 1e4, 4e4),
            "ddG": pplt.Norm("linear", -1e4, 1e4),
            "rfu": pplt.Norm("linear", 0, 1.0, clip=True),
            "drfu": pplt.Norm("linear", -0.5, 0.5, clip=True),
            "mse": None,
            "foldedness": foldedness_cmap,
        }

        levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        norm = pplt.Norm("segmented", levels=levels)
        self.norms["redundancy"] = pplt.DiscreteNorm(levels=levels, norm=norm)

        levels = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0]
        norm = pplt.Norm("segmented", levels=levels)
        self.norms["resolution"] = pplt.DiscreteNorm(levels=levels, norm=norm)

        self.cmaps = {
            "dG": set_bad(pplt.Colormap(tol_cmap("rainbow_PuRd")).reversed()),
            "ddG": set_bad(tol_cmap("PRGn").reversed()),
            "rfu": set_bad(pplt.Colormap(cc.cm.gouldian)),
            "drfu": set_bad(pplt.Colormap(cc.cm.diverging_bwr_20_95_c54)),
            "mse": set_bad(pplt.Colormap("cividis"), color="#e3e3e3"),
            "foldedness": foldedness_cmap,
        }

        colors = ["#6EA72A", "#DAD853", "#FFA842", "#A22D46", "#5D0496"][::-1]
        cmap_redundancy = pplt.Colormap(
            colors, discrete=True, N=len(colors), listmode="discrete"
        )
        cmap_redundancy.set_over("#0E4A21")
        cmap_redundancy.set_bad(NO_COVERAGE)
        self.cmaps["redundancy"] = cmap_redundancy

        colors = ["#008832", "#72D100", "#FFFF04", "#FFB917", "#FF8923"]
        cmap_redundancy = pplt.Colormap(
            colors, discrete=True, N=len(colors), listmode="discrete"
        )
        cmap_redundancy.set_over("#FE2B2E")
        cmap_redundancy.set_bad(NO_COVERAGE)
        self.cmaps["resolution"] = cmap_redundancy

    def __getitem__(self, item):
        cmap = self.cmaps[item]
        norm = self.norms[item]
        return cmap, norm

    def get(self, item, default=None):
        try:
            return self[item]
        except KeyError:
            return default


CMAP_NORM_DEFAULTS = ColorTransforms()


def pymol_figures(
    data,
    output_path,
    pdb_file,
    reference=None,
    field="dG",
    cmap=None,
    norm=None,
    extent=None,
    orient=True,
    views=None,
    name_suffix="",
    additional_views=None,
    img_size=(640, 640),
):

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
        plot_data.columns = pd.MultiIndex.from_product(
            [plot_data.columns, [field]], names=["State", "quantity"]
        )

        cmap_default, norm_default = CMAP_NORM_DEFAULTS["ddG"]
    else:
        plot_data = data
        cmap_default, norm_default = CMAP_NORM_DEFAULTS["dG"]

    cmap = cmap or cmap_default
    norm = norm or norm_default

    for state in protein_states:
        if state == reference_state:
            continue

        values = plot_data[state, field]
        rmin, rmax = extent or [None, None]
        rmin = rmin or values.index.min()
        rmax = rmax or values.index.max()

        values = values.reindex(pd.RangeIndex(rmin, rmax + 1, name="r_number"))
        colors = apply_cmap(values, cmap, norm)
        name = (
            f"pymol_ddG_{state}_ref_{reference_state}"
            if reference_state
            else f"pymol_dG_{state}"
        )
        name += name_suffix
        pymol_render(
            output_path,
            pdb_file,
            colors,
            name=name,
            orient=orient,
            views=views,
            additional_views=additional_views,
            img_size=img_size,
        )


def pymol_render(
    output_path,
    pdb_file,
    colors,
    name="Pymol render",
    orient=True,
    views=None,
    additional_views=None,
    img_size=(640, 640),
):
    if cmd is None:
        raise ModuleNotFoundError("Pymol module is not installed")

    px, py = img_size

    cmd.reinitialize()
    cmd.load(pdb_file)
    if orient:
        cmd.orient()
    cmd.set("antialias", 2)
    cmd.set("fog", 0)

    color_pymol(colors, cmd)

    if views:
        for i, view in enumerate(views):
            cmd.set_view(view)
            cmd.ray(px, py, renderer=0, antialias=2)
            output_file = output_path / f"{name}_view_{i}.png"
            cmd.png(str(output_file))

    else:
        cmd.ray(px, py, renderer=0, antialias=2)
        output_file = output_path / f"{name}_xy.png"
        cmd.png(str(output_file))

        cmd.rotate("x", 90)

        cmd.ray(px, py, renderer=0, antialias=2)
        output_file = output_path / f"{name}_xz.png"
        cmd.png(str(output_file))

        cmd.rotate("z", -90)

        cmd.ray(px, py, renderer=0, antialias=2)
        output_file = output_path / f"{name}_yz.png"
        cmd.png(str(output_file))

        additional_views = additional_views or []

        for i, view in enumerate(additional_views):
            cmd.set_view(view)
            cmd.ray(px, py, renderer=0, antialias=2)
            output_file = output_path / f"{name}_view_{i}.png"
            cmd.png(str(output_file))


def add_cbar(ax, cmap, norm, **kwargs):
    """Truncate or expand cmap such that it covers axes limit and and colorbar to axes"""

    N = cmap.N
    ymin, ymax = np.min(ax.get_ylim()), np.max(ax.get_ylim())
    values = np.linspace(ymin, ymax, num=N)

    norm_clip = copy(norm)
    norm_clip.clip = True
    colors = cmap(norm_clip(values))

    if isinstance(cmap, pplt.DiscreteColormap):
        listmode = "discrete"
    elif isinstance(cmap, pplt.ContinuousColormap):
        listmode = "continuous"
    else:
        listmode = "perceptual"

    cb_cmap = pplt.Colormap(colors, listmode=listmode)

    cb_norm = pplt.Norm("linear", vmin=ymin, vmax=ymax)  # todo allow log norms?
    cbar_kwargs = {**CBAR_KWARGS, **kwargs}
    reverse = np.diff(ax.get_ylim()) < 0

    cbar = ax.colorbar(cb_cmap, norm=cb_norm, reverse=reverse, **cbar_kwargs)

    return cbar


# https://stackoverflow.com/questions/38629830/how-to-turn-off-autoscaling-in-matplotlib-pyplot
@contextmanager
def autoscale_turned_off(ax=None):
    ax = ax or plt.gca()
    lims = [ax.get_xlim(), ax.get_ylim()]
    yield
    ax.set_xlim(*lims[0])
    ax.set_ylim(*lims[1])


def stripplot(
    data,
    ax=None,
    jitter=0.25,
    colors=None,
    offset=0.0,
    orientation="vertical",
    **scatter_kwargs,
):
    ax = ax or plt.gca()
    color_list = _prepare_colors(colors, len(data))

    for i, (d, color) in enumerate(zip(data, color_list)):
        jitter_offsets = (np.random.rand(d.size) - 0.5) * jitter
        cat_var = (
            i * np.ones_like(d) + jitter_offsets + offset
        )  # categorical axis variable
        if orientation == "vertical":
            ax.scatter(cat_var, d, color=color, **scatter_kwargs)
        elif orientation == "horizontal":
            ax.scatter(d, len(data) - cat_var, color=color, **scatter_kwargs)


def _prepare_colors(colors, N):
    if not isinstance(colors, list):
        return [colors] * N
    else:
        return colors


# From joyplot
def _x_range(data, extra=0.2):
    """Compute the x_range, i.e., the values for which the
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
    return np.linspace(
        np.nanmin(data) - extra * sample_range,
        np.nanmax(data) + extra * sample_range,
        1000,
    )


def kdeplot(
    data,
    ax=None,
    offset=0.0,
    orientation="vertical",
    linecolor=None,
    linewidth=None,
    zero_line=True,
    x_extend=1e-3,
    y_scale=None,
    y_norm=None,
    fillcolor=False,
    fill_cmap=None,
    fill_norm=None,
):
    assert not (y_scale and y_norm), "Cannot set both 'y_scale' and 'y_norm'"
    y_scale = 1.0 if y_scale is None else y_scale

    color_list = _prepare_colors(linecolor, len(data))

    for i, (d, color) in enumerate(zip(data, color_list)):
        # todo remove NaNs?

        # Perhaps also borrow this part from joyplot
        kde_func = kde.gaussian_kde(d)
        kde_x = _x_range(d, extra=0.4)
        kde_y = kde_func(kde_x) * y_scale
        if y_norm:
            kde_y = y_norm * kde_y / kde_y.max()
        bools = kde_y > x_extend * kde_y.max()
        kde_x = kde_x[bools]
        kde_y = kde_y[bools]

        cat_var = len(data) - i + kde_y + offset  # x in horizontal
        cat_var_zero = (len(data) - i) * np.ones_like(kde_y) + offset

        # x = i * np.ones_like(d) + jitter_offsets + offset  # 'x' like, could be y axis
        if orientation == "horizontal":
            plot_x = kde_x
            plot_y = cat_var
            img_data = kde_x.reshape(1, -1)
        elif orientation == "vertical":
            plot_x = len(data) - cat_var
            plot_y = kde_x
            img_data = kde_x[::-1].reshape(-1, 1)
        else:
            raise ValueError(f"Invalid value '{orientation}' for 'orientation'")

        (line,) = ax.plot(plot_x, plot_y, color=color, linewidth=linewidth)
        if zero_line:
            ax.plot(
                [plot_x[0], plot_x[-1]],
                [plot_y[0], plot_y[-1]],
                color=line.get_color(),
                linewidth=linewidth,
            )

        if fillcolor:
            # todo refactor to one if/else orientation
            color = line.get_color() if fillcolor is True else fillcolor
            if orientation == "horizontal":
                ax.fill_between(
                    kde_x,
                    plot_y,
                    np.linspace(plot_y[0], plot_y[-1], num=plot_y.size, endpoint=True),
                    color=color,
                )
            elif orientation == "vertical":
                ax.fill_betweenx(
                    kde_x, len(data) - cat_var, len(data) - cat_var_zero, color=color
                )

        if fill_cmap:
            fill_norm = fill_norm or pplt.Norm("linear")

            xmin, xmax = np.min(plot_x), np.max(plot_x)
            ymin, ymax = np.min(plot_y), np.max(plot_y)
            extent = (
                [xmin - offset, xmax - offset, ymin, ymax]
                if orientation == "horizontal"
                else [xmin, xmax, ymin - offset, ymax - offset]
            )
            im = Axes.imshow(
                ax,
                img_data,
                aspect="auto",
                cmap=fill_cmap,
                norm=fill_norm,
                extent=extent,
            )  # left, right, bottom, top
            (fill_line,) = ax.fill(plot_x, plot_y, facecolor="none")
            im.set_clip_path(fill_line)


def boxplot(
    data,
    ax,
    offset=0.0,
    orientation="vertical",
    widths=0.25,
    linewidth=None,
    linecolor=None,
    **kwargs,
):
    if orientation == "vertical":
        vert = True
        positions = np.arange(len(data)) + offset
    elif orientation == "horizontal":
        vert = False
        positions = len(data) - np.arange(len(data)) - offset
    else:
        raise ValueError(
            f"Invalid value '{orientation}' for 'orientation', options are 'horizontal' or 'vertical'"
        )

    # todo for loop
    boxprops = kwargs.pop("boxprops", {})
    whiskerprops = kwargs.pop("whiskerprops", {})
    medianprops = kwargs.pop("whiskerprops", {})

    boxprops["linewidth"] = linewidth
    whiskerprops["linewidth"] = linewidth
    medianprops["linewidth"] = linewidth

    boxprops["color"] = linecolor
    whiskerprops["color"] = linecolor
    medianprops["color"] = linecolor

    Axes.boxplot(
        ax,
        data,
        vert=vert,
        positions=positions,
        widths=widths,
        boxprops=boxprops,
        whiskerprops=whiskerprops,
        medianprops=medianprops,
        **kwargs,
    )


def label_axes(labels, ax, offset=0.0, orientation="vertical", **kwargs):
    # todo check offset sign
    if orientation == "vertical":
        ax.set_xticks(np.arange(len(labels)) + offset)
        ax.set_xticklabels(labels, **kwargs)
    elif orientation == "horizontal":
        ax.set_yticks(len(labels) - np.arange(len(labels)) + offset)
        ax.set_yticklabels(labels, **kwargs)


class FitResultPlotBase(object):
    def __init__(self, fit_result):
        self.fit_result = fit_result

    # todo equivalent this for axes?
    def _make_figure(self, figure_name, **kwargs):
        if not figure_name.endswith("_figure"):
            figure_name += "_figure"

        function = globals()[figure_name]
        args_dict = self._get_arg(figure_name)

        # return dictionary
        # keys: either protein state name (hdxm.name) or 'All states'

        figures_dict = {
            name: function(arg, **kwargs) for name, arg in args_dict.items()
        }
        return figures_dict

    def make_figure(self, figure_name, **kwargs):
        figures_dict = self._make_figure(figure_name, **kwargs)
        if len(figures_dict) == 1:
            return next(iter(figures_dict.values()))
        else:
            return figures_dict

    def get_fit_timepoints(self):
        all_timepoints = np.concatenate(
            [hdxm.timepoints for hdxm in self.fit_result.hdxm_set]
        )

        # x_axis_type = self.settings.get('fit_time_axis', 'Log')
        x_axis_type = "Log"  # todo configureable
        num = 100
        if x_axis_type == "Linear":
            time = np.linspace(0, all_timepoints.max(), num=num)
        elif x_axis_type == "Log":
            elem = all_timepoints[np.nonzero(all_timepoints)]
            start = np.log10(elem.min())
            end = np.log10(elem.max())
            pad = (end - start) * 0.1
            time = np.logspace(start - pad, end + pad, num=num, endpoint=True)
        else:
            raise ValueError("Invalid value for 'x_axis_type'")

        return time

    # repeated code with fitreport (pdf) -> base class for fitreport
    def _get_arg(self, plot_func_name):
        # Add _figure suffix if not present
        if not plot_func_name.endswith("_figure"):
            plot_func_name += "_figure"

        if plot_func_name == "peptide_coverage_figure":
            return {hdxm.name: hdxm.data for hdxm in self.fit_result.hdxm_set.hdxm_list}
        elif plot_func_name == "residue_time_scatter_figure":
            return {hdxm.name: hdxm for hdxm in self.fit_result.hdxm_set.hdxm_list}
        elif plot_func_name == "residue_scatter_figure":
            return {"All states": self.fit_result.hdxm_set}
        elif plot_func_name == "dG_scatter_figure":
            return {"All states": self.fit_result.output}
        elif plot_func_name == "ddG_scatter_figure":
            return {"All states": self.fit_result.output}
        elif plot_func_name == "linear_bars_figure":
            return {"All states": self.fit_result.output}
        elif plot_func_name == "rainbowclouds_figure":
            return {"All states": self.fit_result.output}
        elif plot_func_name == "peptide_mse_figure":
            return {"All states": self.fit_result.get_peptide_mse()}
        elif plot_func_name == "loss_figure":
            return {"All states": self.fit_result}
        else:
            raise ValueError(f"Unknown plot function {plot_func_name!r}")


ALL_PLOT_TYPES = [
    "peptide_coverage",
    "residue_scatter",
    "dG_scatter",
    "ddG_scatter",
    "linear_bars",
    "rainbowclouds",
    "peptide_mse",
    "loss",
]


class FitResultPlot(FitResultPlotBase):
    def __init__(self, fit_result, output_path=None, **kwargs):
        super().__init__(fit_result)
        self.output_path = Path(output_path) if output_path else None
        self.output_path.mkdir(exist_ok=True)
        if self.output_path and not self.output_path.is_dir():
            raise ValueError(f"Output path {output_path!r} is not a valid directory")

        # todo save kwargs / rc params? / style context (https://matplotlib.org/devdocs/tutorials/introductory/customizing.html)

    def save_figure(self, fig_name, ext=".png", **kwargs):
        figures_dict = self._make_figure(fig_name, **kwargs)

        if self.output_path is None:
            raise ValueError(
                f"No output path given when `FitResultPlot` object as initialized"
            )
        for name, fig_tup in figures_dict.items():
            fig = fig_tup if isinstance(fig_tup, plt.Figure) else fig_tup[0]

            if name == "All states":  # todo variable for 'All states'
                file_name = f"{fig_name.replace('_figure', '')}{ext}"
            else:
                file_name = f"{fig_name.replace('_figure', '')}_{name}{ext}"
            file_path = self.output_path / file_name
            fig.savefig(file_path)
            plt.close(fig)

    def plot_all(self, **kwargs):
        pbar = tqdm(ALL_PLOT_TYPES)
        for plot_type in pbar:
            pbar.set_description(plot_type)
            fig_kwargs = kwargs.get(plot_type, {})
            self.save_figure(plot_type, **fig_kwargs)


def plot_fitresults(
    fitresult_path,
    reference=None,
    plots="all",
    renew=False,
    cmap_and_norm=None,
    output_path=None,
    output_type=".png",
    **save_kwargs,
):
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

    raise DeprecationWarning(
        "This function is deprecated, use FitResultPlot.plot_all instead"
    )
    # batch results only
    history_path = fitresult_path / "model_history.csv"
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

    # todo needs tidying up
    cmap_and_norm = cmap_and_norm or {}
    dG_cmap, dG_norm = cmap_and_norm.get("dG", (None, None))
    dG_cmap_default, dG_norm_default = default_cmap_norm("dG")
    ddG_cmap, ddG_norm = cmap_and_norm.get("ddG", (None, None))
    ddG_cmap_default, ddG_norm_default = default_cmap_norm("ddG")
    dG_cmap = ddG_cmap or dG_cmap_default
    dG_norm = dG_norm or dG_norm_default
    ddG_cmap = ddG_cmap or ddG_cmap_default
    ddG_norm = ddG_norm or ddG_norm_default

    # check_exists = lambda x: False if renew else x.exists()
    # todo add logic for checking renew or not

    if plots == "all":
        plots = [
            "loss",
            "rfu_coverage",
            "rfu_scatter",
            "dG_scatter",
            "ddG_scatter",
            "linear_bars",
            "rainbowclouds",
            "peptide_mse",
        ]

    # def check_update(pth, fname, extensions, renew):
    #     # Returns True if the target graph should be renewed or not
    #     if renew:
    #         return True
    #     else:
    #         pths = [pth / (fname + ext) for ext in extensions]
    #         return any([not pth.exists() for pth in pths])

    # plots = [p for p in plots if check_update(output_path, p, output_type, renew)]

    if "loss" in plots:
        loss_df = fitresult.losses
        loss_df.plot()

        mse_loss = loss_df["mse_loss"]
        reg_loss = loss_df.iloc[:, 1:].sum(axis=1)
        reg_percentage = 100 * reg_loss / (mse_loss + reg_loss)
        fig = plt.gcf()
        ax = plt.gca()
        ax1 = ax.twinx()
        reg_percentage.plot(ax=ax1, color="k")
        ax1.set_xlim(0, None)
        for ext in output_type:
            f_out = output_path / ("loss" + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if "rfu_coverage" in plots:
        for hdxm in fitresult.hdxm_set:
            fig, axes, cbar_ax = peptide_coverage_figure(hdxm.data)
            for ext in output_type:
                f_out = output_path / (f"rfu_coverage_{hdxm.name}" + ext)
                plt.savefig(f_out)
            plt.close(fig)

    # todo rfu_scatter_timepoint

    if "rfu_scatter" in plots:
        fig, axes, cbar = residue_scatter_figure(fitresult.hdxm_set)
        for ext in output_type:
            f_out = output_path / (f"rfu_scatter" + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if "dG_scatter" in plots:
        fig, axes, cbars = dG_scatter_figure(
            fitresult.output.df, cmap=dG_cmap, norm=dG_norm
        )
        for ext in output_type:
            f_out = output_path / (f"dG_scatter" + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if "ddG_scatter" in plots:
        fig, axes, cbars = ddG_scatter_figure(
            fitresult.output.df, reference=reference, cmap=ddG_cmap, norm=ddG_norm
        )
        for ext in output_type:
            f_out = output_path / (f"ddG_scatter" + ext)
            plt.savefig(f_out)
        plt.close(fig)

    if "linear_bars" in plots:
        fig, axes = linear_bars_figure(fitresult.output.df)
        for ext in output_type:
            f_out = output_path / (f"dG_linear_bars" + ext)
            plt.savefig(f_out)
        plt.close(fig)

        if reference_state:
            fig, axes = linear_bars_figure(fitresult.output.df, reference=reference)
            for ext in output_type:
                f_out = output_path / (f"ddG_linear_bars" + ext)
                plt.savefig(f_out)
            plt.close(fig)

    if "rainbowclouds" in plots:
        fig, ax = rainbowclouds_figure(fitresult.output.df)
        for ext in output_type:
            f_out = output_path / (f"dG_rainbowclouds" + ext)
            plt.savefig(f_out)
        plt.close(fig)

        if reference_state:
            fig, axes = rainbowclouds_figure(fitresult.output.df, reference=reference)
            for ext in output_type:
                f_out = output_path / (f"ddG_rainbowclouds" + ext)
                plt.savefig(f_out)
            plt.close(fig)

    if "peptide_mse" in plots:
        fig, axes, cbars = peptide_mse_figure(fitresult.get_peptide_mse())
        for ext in output_type:
            f_out = output_path / (f"peptide_mse" + ext)
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
