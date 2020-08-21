import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import uuid
import shutil
from functools import lru_cache, partial
from pyhdx.support import grouper, autowrap
from pyhdx.plot import plot_peptides
from tqdm.auto import tqdm

try:
    import pylatex as pyl

except ModuleNotFoundError:
    pyl = None
import tempfile


geometry_options = {
    "lmargin": "1in",
    "rmargin": "1.5in"
    }


class Report(object):
    """

    .pdf output document
    """

    def __init__(self, output, name=None, doc=None, add_date=True):
        if not pyl:
            raise ModuleNotFoundError('pylatex module not installed')
        name = name or output.series.state
        self.output = output
        self.doc = doc or self._init_doc(name, add_date=add_date)
        self._temp_dir = self.make_temp_dir()

    def make_temp_dir(self):
        _tmp_path = os.path.abspath(os.path.join(tempfile.gettempdir(), str(id(self))))

        if not os.path.exists(_tmp_path):
            os.makedirs(_tmp_path)
        return _tmp_path

    def _init_doc(self, name, add_date=True):
        doc = pyl.Document(name, geometry_options=geometry_options)
        doc.packages.append(pyl.Package('hyperref'))
        doc.preamble.append(pyl.Command('title', f'Supplementary Figures for {name}'))
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
        funcs = [partial(self.output._make_coverage_graph, i, **kwargs) for i in range(len(self.output.series))]
        self.make_subfigure(funcs, layout=layout, close=close)

    def add_peptide_figures(self, layout=(5, 4), close=True, **kwargs):
        funcs = [partial(self.output._make_peptide_graph, i, **kwargs) for i in range(len(self.output.series.cov))]
        self.make_subfigure(funcs, layout=layout, close=close)

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

    def generate_pdf(self, file_path=None):
        self.doc.generate_pdf(file_path, compiler='pdflatex')


class Output(object):

    def __init__(self, series, fit_results, **settings):
        #todo series is now an attribute in fit_results
        self.series = series
        self.fit_results = fit_results # dictionary with name: KineticsFitresult
        self.settings = settings #dictionary with additional settings

    @property
    def fit_timepoints(self):
        timepoints = self.series.timepoints
        x_axis_type = self.settings.get('fit_time_axis', 'Log')
        if x_axis_type == 'Linear':
            time = np.linspace(0, timepoints.max(), num=250)
        elif x_axis_type == 'Log':
            time = np.logspace(-2, np.log10(timepoints.max()), num=250)
        return time

    @lru_cache(8)
    def call_fitresult(self, fit_name, time=None):
        """

        Parameters
        ----------
        fit_name: :obj:`str`
            Name of the fit result to use from `fit_results` dictionary
        time: array_like
            Optional array of timepoints

        Returns
        -------
        uptake: :obj:`str`
            Numpy array with fitted uptake values for each peptide per row
        """
        fit_timepoints = time or self.fit_timepoints

        # Moved to fit results
        # fit_result = self.fit_results[fit_name]
        # d_list = []
        # if fit_result.model_type == 'Single':
        #     for time in fit_timepoints:
        #         p = fit_result.get_p(time)
        #         p = np.nan_to_num(p)
        #         d = self.series.cov.X.dot(p)
        #         d_list.append(d)
        # elif fit_result.model_type == 'Global':
        #     for time in fit_timepoints:
        #         d = fit_result.get_d(time)
        #         d_list.append(d)
        #
        # uptake = np.vstack(d_list).T

        uptake = self.fit_results[fit_name](fit_timepoints)

        return uptake

    def add_peptide_fits(self, ax_scale='log', fit_names=None):
        pass

    def peptide_graph_generator(self, **kwargs):
        for i in range(len(self.series.cov)):
            yield from self._make_peptide_graph(i, **kwargs)

    def _make_peptide_graph(self, index, figsize=(4,4), ax_scale='log', **fig_kwargs):
        """yield single peptide grpahs"""

        fig, ax = plt.subplots(figsize=figsize)
        if ax_scale == 'log':
            ax.set_xscale('log')

        for k in self.fit_results.keys():
            ax.plot(self.fit_timepoints, self.call_fitresult(k)[index], label=k)

        ax.scatter(self.series.timepoints, self.series.uptake_corrected.T[index], color='k')
        t_unit = fig_kwargs.get('time_unit', '')
        t_unit = f'({t_unit})' if t_unit else t_unit
        ax.set_xlabel(f'Time' + t_unit)
        ax.set_ylabel('Corrected D-uptake')
        start, end = self.series.cov.data['_start'][index], self.series.cov.data['_end'][index]
        ax.set_title(f'peptide_{start}_{end}')
        ax.legend()
        plt.tight_layout()

        return fig

    def _make_coverage_graph(self, index, figsize=(14, 4), cbar=True, **fig_kwargs):
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
