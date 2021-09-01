"""
This module allows users to generate a .pdf output report from their HDX measurement

(Currently partially out of date)
"""


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import uuid
import shutil
from functools import lru_cache, partial
from pyhdx.support import grouper, autowrap
from pyhdx.plot import plot_peptides
from pyhdx.fitting_torch import TorchSingleFitResult
from tqdm.auto import tqdm

import pylatex as pyl
import proplot as pplt
import tempfile


geometry_options = {
    "lmargin": "1in",
    "rmargin": "1.5in"
    }


# plot_defaults = {
#     ''}


class Report(object):
    """

    .pdf output document
    """

    def __init__(self, output, title=None, doc=None, add_date=True):
        self.title = title or f'Fit report for {output.fit_result.data_obj.name}'
        self.output = output
        self.doc = doc or self._init_doc(add_date=add_date)
        self._temp_dir = self.make_temp_dir()

    def make_temp_dir(self):
        _tmp_path = os.path.abspath(os.path.join(tempfile.gettempdir(), str(id(self))))

        if not os.path.exists(_tmp_path):
            os.makedirs(_tmp_path)
        return _tmp_path

    def _init_doc(self, add_date=True):
        doc = pyl.Document(geometry_options=geometry_options)
        doc.packages.append(pyl.Package('hyperref'))
        doc.preamble.append(pyl.Command('title', self.title))
        if add_date:
            doc.preamble.append(pyl.Command('date', pyl.NoEscape(r'\today')))
        else:
            doc.preamble.append(pyl.Command('date', pyl.NoEscape(r'')))
        doc.append(pyl.NoEscape(r'\maketitle'))
        doc.append(pyl.NewPage())

        return doc

    def _save_fig(self, fig, *args, extension='pdf', **kwargs):
        filename = '{}.{}'.format(str(uuid.uuid4()), extension.strip('.'))
        filepath = os.path.join(self._temp_dir, filename)
        fig.savefig(filepath, *args, **kwargs)
        return filepath

    def test_mpl(self):
        fig = plt.figure()
        plt.plot([2,3,42,1])

        file_path = self._save_fig(fig)

        with self.doc.create(pyl.Figure(position='htbp')) as plot:
            plot.add_image(pyl.NoEscape(file_path), width=pyl.NoEscape(r'1\textwidth'))
            plot.add_caption('I am a caption.')

    def add_coverage_figures(self, layout=(6, 2), close=True, **kwargs):
        raise NotImplementedError()
        funcs = [partial(self.output._make_coverage_graph, i, **kwargs) for i in range(len(self.output.series))]
        self.make_subfigure(funcs, layout=layout, close=close)

    def add_peptide_figures(self, ncols=4, nrows=5, **kwargs):

        Np = self.output.fit_result.data_obj.Np
        indices = range(Np)
        n = ncols*nrows
        chunks = [indices[i:i + n] for i in range(0, len(indices), n)]
        for chunk in tqdm(chunks):
            with self.doc.create(pyl.Figure(position='ht')) as tex_fig:
                fig = self.output._make_peptide_subplots(chunk, ncols=ncols, nrows=nrows, **kwargs)
                file_path = self._save_fig(fig)
                plt.close(fig)

                tex_fig.add_image(file_path, width=pyl.NoEscape(r'\textwidth'))

        #self.make_subfigure(funcs, layout=layout, close=close)

    def make_subfigure(self, fig_funcs, layout=(5, 4), close=True):
        #todo figure out how to iterate properly
        n = np.product(layout)
        chunks = grouper(n, fig_funcs)
        w = str(1/layout[1])
        pbar = tqdm(total=len(fig_funcs))
        for chunk in chunks:
            with self.doc.create(pyl.Figure(position='ht')) as tex_fig:
                for i, fig_func in enumerate(chunk):
                    if fig_func is None:
                        continue
                    with self.doc.create(pyl.SubFigure(position='b', width=pyl.NoEscape(w + r'\linewidth'))) as subfig:
                        fig = fig_func()
                        file_path = self._save_fig(fig, bbox_inches='tight') # todo access these kwargs
                        if close:
                            plt.close(fig)
                        subfig.add_image(file_path, width=pyl.NoEscape(r'\linewidth'))
                    if i % layout[1] == layout[1] - 1:
                        self.doc.append('\n')
                    pbar.update(1)

            self.doc.append(pyl.NewPage())

    def test_subfigure(self):
        fig = plt.figure()
        plt.plot([2,3,42,1])

        file_path = self._save_fig(fig)

        with self.doc.create(pyl.Figure(position='h!')) as kittens:
            w = str(0.25)
            for i in range(8):
                with self.doc.create(pyl.SubFigure(
                        position='b',
                        width=pyl.NoEscape(w + r'\linewidth'))) as left_kitten:

                    left_kitten.add_image(file_path,
                                          width=pyl.NoEscape(r'\linewidth'))
                    left_kitten.add_caption(f'Kitten on the {i}')
                if i % 4 == 3:
                    self.doc.append('\n')
            kittens.add_caption("Two kittens")

    def rm_temp_dir(self):
        """Remove the temporary directory specified in ``_tmp_path``."""

        if os.path.exists(self._temp_dir):
            shutil.rmtree(self._temp_dir)

    def generate_pdf(self, file_path):
        self.doc.generate_pdf(file_path, compiler='pdflatex')


class Output(object):
    # Currently only TorchSingleFitResult support
    def __init__(self, fit_result, time_axis=None, **settings):
        assert isinstance(fit_result, TorchSingleFitResult), "Invalid type of `fit_result`"
        self.settings = {'fit_time_axis': 'Log'}
        self.settings.update(settings)

        #todo restore multiple fit results functionality
        self.fit_result = fit_result
        self.fit_timepoints = time_axis or self.get_fit_timepoints()
        self.d_calc = self.fit_result(self.fit_timepoints)

    def get_fit_timepoints(self):
        timepoints = self.fit_result.data_obj.timepoints
        x_axis_type = self.settings.get('fit_time_axis', 'Log')
        num = 100
        if x_axis_type == 'Linear':
            time = np.linspace(0, timepoints.max(), num=num)
        elif x_axis_type == 'Log':
            elem = timepoints[np.nonzero(timepoints)]
            time = np.logspace(np.log10(elem.min()) - 1, np.log10(elem.max()), num=num, endpoint=True)

        return time

    def add_peptide_fits(self, ax_scale='log', fit_names=None):
        pass

    def peptide_graph_generator(self, **kwargs):
        for i in range(len(self.series.coverage)):
            yield from self._make_peptide_graph(i, **kwargs)

    def _make_peptide_subplots(self, indices, **fig_kwargs):
        """yield single peptide grpahs"""
        nrows = fig_kwargs.pop('nrows', int(np.floor(np.sqrt(len(indices)))))
        ncols = fig_kwargs.pop('ncols', int(np.ceil(len(indices) / nrows)))

        default_kwargs = {'sharex': 1, 'sharey': 1, 'ncols': ncols, 'nrows': nrows}
        default_kwargs.update(fig_kwargs)

        fig, axes = pplt.subplots(**default_kwargs)
        axes_iter = iter(axes)
        for i, ax in zip(indices, axes_iter):
            ax.plot(self.fit_timepoints, self.d_calc[i], color='r')
            ax.scatter(self.fit_result.data_obj.timepoints, self.fit_result.data_obj.d_exp[i], color='k')

            start = self.fit_result.data_obj.coverage.data['_start'][i]
            end = self.fit_result.data_obj.coverage.data['_end'][i]
            ax.set_title(f'Peptide_{i}: {start} - {end}')

        t_unit = fig_kwargs.get('time_unit', 'min')
        t_unit = f'({t_unit})' if t_unit else t_unit

        # turn off remaining axes
        #todo proplot issue
        axes.format(xscale='log', xlabel=f'Time' + t_unit, ylabel='Corrected D-uptake', xformatter='log')
        xlim = axes[0].get_xlim()
        for ax in axes_iter:
            #ax.axis('off')
            ax.set_axis_off()
        axes.format(xlim=xlim)

        return fig

    def _make_peptide_graph(self, index, figsize=(4,4), ax_scale='log', **fig_kwargs):
        """yield single peptide grpahs"""

        fig, ax = plt.subplots(figsize=figsize)
        if ax_scale == 'log':
            ax.set_xscale('log')
            ax.get_xaxis().get_major_formatter().set_scientific(True)

        ax.plot(self.fit_timepoints, self.d_calc[index], color='r')
        ax.scatter(self.fit_result.data_obj.timepoints, self.fit_result.data_obj.d_exp[index], color='k')

        t_unit = fig_kwargs.get('time_unit', 'min')
        t_unit = f'({t_unit})' if t_unit else t_unit
        ax.set_xlabel(f'Time' + t_unit)
        ax.set_ylabel('Corrected D-uptake')
        start = self.fit_result.data_obj.coverage.data['_start'][index]
        end = self.fit_result.data_obj.coverage.data['_end'][index]
        ax.set_title(f'peptide_{start}_{end}')


        #ax.legend()
        plt.tight_layout()

        return fig

    def _make_coverage_graph(self, index, figsize=(14, 4), cbar=True, **fig_kwargs):
        raise NotImplementedError("coverage not implemented")
        peptides = self.series[index]

        cmap = fig_kwargs.get('cmap', 'jet')
        if cbar:
            fig, (ax_main, ax_cbar) = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [40, 1], 'wspace': 0.025})

            norm = mpl.colors.Normalize(vmin=0, vmax=100)
            cmap = mpl.cm.get_cmap(cmap)
            cb1 = mpl.colorbar.ColorbarBase(ax_cbar, cmap=mpl.cm.get_cmap(cmap),
                                            norm=norm,
                                            orientation='vertical', ticks=[0, 100])
            cb1.set_label('Uptake %', x=-1, rotation=270)
            # cbar_ax.xaxis.set_ticks_position('top')
            cb1.set_ticks([0, 100])

        else:
            fig, ax_main = plt.subplots(figsize=figsize)

        wrap = autowrap(peptides)
        plot_peptides(peptides, wrap, ax_main, **fig_kwargs)
        ax_main.set_xlabel('Residue number')
        t_unit = fig_kwargs.get('time_unit', '')
        fig.suptitle(f'Deuterium uptake at t={peptides.exposure} ' + t_unit)
        plt.tight_layout()

        return fig
