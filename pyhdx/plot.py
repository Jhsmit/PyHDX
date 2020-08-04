import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from bokeh.models import LinearColorMapper, ColorBar, ColumnDataSource, Rect, LabelSet, HoverTool
from bokeh.plotting import figure
from bokeh.layouts import row, column


def plot_residue_map(pm, scores=None, ax=None, cmap='jet', bad='k', cbar=True, **kwargs):
    img = (pm.X > 0).astype(float)
    if scores is not None:
        img *= scores[:, np.newaxis]
    elif pm.scores is not None:
        img *= pm.scores[:, np.newaxis]

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


def plot_peptides(pm, wrap, ax, color=True, labels=False, cmap='jet', **kwargs):
    """
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

    rect_kwargs = {'linewidth': 1, 'linestyle': '-', 'edgecolor': 'k'}
    rect_kwargs.update(kwargs)

    cmap = mpl.cm.get_cmap(cmap)
    i = -1

    for p_num, e in enumerate(pm.data):
        if i < -wrap:
            i = -1

        if color:
            c = cmap(e['scores'] / 100)
        else:
            c = '#707070'

        width = e['end'] - e['start'] + 1
        rect = Rectangle((e['start'] - 0.5, i), width, 1, facecolor=c, **rect_kwargs)
        ax.add_patch(rect)
        if labels:
            rx, ry = rect.get_xy()
            cy = ry# + rect.get_height() / 2.
            cx = rx
            ax.annotate(str(p_num), (cx, cy), color='k', fontsize=6, va='bottom', ha='right')

        i -= 1

    ax.set_ylim(-wrap, 0)
    ax.set_xlim(0, pm.end + 5)
    ax.set_yticks([])


def make_coverage_figure(pm, wrap, aa_per_subplot, color=False, figsize=(10, 8), labels=False, **kwargs):
    rect_kwargs = {'linewidth': 1, 'linestyle': '-', 'edgecolor': 'k'}
    rect_kwargs.update(kwargs)

    num_axes = pm.end // aa_per_subplot + 1
    cmap = mpl.cm.get_cmap('viridis')

    fig, axes = plt.subplots(num_axes, figsize=figsize)
    axes = [axes] if num_axes == 1 else axes
    for j, ax in enumerate(axes):
        i = -1

        for p_num, e in enumerate(pm.data):
            if i < -wrap:
                i = -1

            if color:
                c = cmap(e['scores'] / 100)
            else:
                c = '#313695'

            width = e['end'] - e['start'] + 1
            rect = Rectangle((e['start'] - 0.5, i), width, 1, facecolor=c, **rect_kwargs)

            ax.add_patch(rect)
            if labels:
                rx, ry = rect.get_xy()
                cy = ry# + rect.get_height() / 2.
                cx = rx
                ax.annotate(str(p_num), (cx, cy), color='k', fontsize=6, va='bottom', ha='right')

            i -= 1

        ax.set_ylim(-wrap, 0)
        ax.set_xlim(j * aa_per_subplot, (j + 1) * aa_per_subplot)
        ax.set_yticks([])

    if color:
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.87, 0.05, 0.02, 0.9])
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=mpl.cm.get_cmap(cmap),
                                    norm=norm,
                                    orientation='vertical')
        cb1.set_label('Uptake (%)')
    else:
        plt.tight_layout()

    axes[-1].set_xlabel('Residue number')


    return fig, axes


def _bokeh_coverage(pm, wrap, aa_per_subplot, color=False, labels=False, **kwargs):
    num_axes = pm.end // aa_per_subplot + 1
    cmap = mpl.cm.get_cmap('jet')
    c_rgba = cmap(pm.data['scores'] / 100)
    c = [mpl.colors.to_hex(color) for color in c_rgba]

    pal = tuple(mpl.colors.to_hex(cmap(value)) for value in np.linspace(0, 1, 256, endpoint=True))
    color_mapper = LinearColorMapper(palette=pal, low=0, high=100)
    color_bar = ColorBar(color_mapper=color_mapper)

    repeats = (len(pm) // wrap) + 1
    y = (list(range(wrap, 0, -1))*repeats)[:len(pm.data)]
    width = pm.data['end'] - pm.data['start'] + 1
    x = pm.data['start'] - 0.5 + (width / 2)
    label_x = pm.data['start']
    names = [str(i) for i in range(len(pm.data))]
    plot_dict = dict(x=x, label_x=label_x, y=y, width=width, c=c, names=names)
    prop_dict = {name: pm.data[name] for name in pm.data.dtype.names}
    source = ColumnDataSource({**plot_dict, **prop_dict})
    glyph = Rect(x='x', y='y', width='width', height=1, fill_color='c')
    labels = LabelSet(x='label_x', y='y', text='names', source=source, text_baseline='middle', text_align='left')

    figures = []
    # plot_width=750, plot_height=int(TOTAL_HEIGHT/num_axes),
    for j in range(num_axes):
        fig = figure(title=None, min_border=0, tools='pan,wheel_zoom,box_zoom,save,reset,hover',
                     x_range=(j*aa_per_subplot, (j + 1) * aa_per_subplot))
        fig.add_glyph(source, glyph)
        fig.add_layout(labels)
        hover = fig.select(dict(type=HoverTool))
        hover.tooltips = [('Pos', '$x{int}'),
                          ('Index', '@names'),
                          ('Start', '@start (@_start)'),
                          ('End', '@end (@_end)'),
                          ('Sequence', '@sequence'),
                          ('Score', '@scores'),
                          ('Uptake', '@uptake (@uptake_corrected / @ex_residues, @maxuptake)')]

        figures.append(fig)

    #https://github.com/bokeh/bokeh/issues/7093
    #ummy = figure(width=100, toolbar_location=None, min_border=0, outline_line_color=None)
    #dummy.add_layout(color_bar, 'right')
    layout = row(column(*figures, sizing_mode='stretch_both'), sizing_mode='stretch_both')

    return layout, figures, labels


"""
coverage = np.sum(s.big_X  >0,axis=0)
tiled = np.tile(coverage, (100, 1))

cmap = mpl.cm.get_cmap('Blues')
cmap.set_under('grey')

def truncate_colormap(cmap, minval=0.5, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


new_cmap = truncate_colormap(cmap, 0.5, 1)
new_cmap.set_under('grey')
norm = mpl.colors.Normalize(vmin=0.25, vmax = tiled.max())
fig, ax = plt.subplots(figsize=(8, 2))
ax.imshow(tiled, cmap=new_cmap, norm=norm)
ax.set_yticks([])
ax.set_xlabel('Residue number')
plt.savefig('Coverage_bar.png', dpi=600)
"""
