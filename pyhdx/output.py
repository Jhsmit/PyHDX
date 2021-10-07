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
from pyhdx.fitting_torch import TorchSingleFitResult
from tqdm.auto import tqdm
from pathlib import Path

import pylatex as pyl
import proplot as pplt
import tempfile

from concurrent import futures

geometry_options = {
    "lmargin": "1in",
    "rmargin": "1.5in"
    }


class BaseReport(object):
    pass


class Report(BaseReport):
    def __init__(self, hdxm_set, **kwargs):
        raise NotImplementedError()


class FitReport(object):
    """
    Create .pdf output of a fit result
    """
    def __init__(self, fit_result, title=None, doc=None, add_date=True, temp_dir=None):
        self.title = title or f'Fit report'
        self.fit_result = fit_result
        self.doc = doc or self._init_doc(add_date=add_date)
        self._temp_dir = temp_dir or self.make_temp_dir()
        self._temp_dir = Path(self._temp_dir)

        self.figure_queue = []
        self.tex_dict = {}  # dictionary gathering lists of partial functions which when executed generate the tex output
        self._figure_number = 0  #todo automate

    def make_temp_dir(self):
        #todo pathlib
        _tmp_path = os.path.abspath(os.path.join(tempfile.gettempdir(), str(id(self))))

        if not os.path.exists(_tmp_path):
            os.makedirs(_tmp_path)
        return _tmp_path

    def _init_doc(self, add_date=True):
        doc = pyl.Document(geometry_options=geometry_options)
        doc.packages.append(pyl.Package('float'))
        doc.packages.append(pyl.Package('hyperref'))

        doc.preamble.append(pyl.Command('title', self.title))
        if add_date:
            doc.preamble.append(pyl.Command('date', pyl.NoEscape(r'\today')))
        else:
            doc.preamble.append(pyl.Command('date', pyl.NoEscape(r'')))
        doc.append(pyl.NoEscape(r'\maketitle'))
        doc.append(pyl.NewPage())
        doc.append(pyl.Command('tableofcontents'))
        doc.append(pyl.NewPage())

        return doc

    def _save_fig(self, fig, *args, extension='pdf', **kwargs):
        filename = '{}.{}'.format(str(uuid.uuid4()), extension.strip('.'))
        filepath = os.path.join(self._temp_dir, filename)
        fig.savefig(filepath, *args, **kwargs)
        return filepath

    def reset_doc(self, add_date=True):
        self.doc = self._init_doc(add_date=add_date)

    def get_fit_timepoints(self):
        all_timepoints = np.concatenate([hdxm.timepoints for hdxm in self.fit_result.data_obj])

        #x_axis_type = self.settings.get('fit_time_axis', 'Log')
        x_axis_type = 'Log' # todo configureable
        num = 100
        if x_axis_type == 'Linear':
            time = np.linspace(0, all_timepoints.max(), num=num)
        elif x_axis_type == 'Log':
            elem = all_timepoints[np.nonzero(all_timepoints)]
            start = np.log10(elem.min())
            end = np.log10(elem.max())
            pad = (end - start)*0.1
            time = np.logspace(start-pad, end+pad, num=num, endpoint=True)
        else:
            raise ValueError("Invalid value for 'x_axis_type'")

        return time

    def figure_number(self):
        self._figure_number += 1
        return self._figure_number

    def add_peptide_uptake_curves(self, layout=(5, 4), time_axis=None):
        extension = '.pdf'
        self.tex_dict['peptide_uptake'] = {}

        nrows, ncols = layout
        n = nrows*ncols
        time = time_axis or self.get_fit_timepoints()
        if time.ndim == 1:
            time = np.tile(time, (len(self.fit_result), 1))

        d_calc = self.fit_result(time)  # Ns x Np x Nt

        fig_factory = partial(pplt.subplots, ncols=ncols, nrows=nrows, sharex=1, sharey=1, num=self.figure_number())

        # iterate over samples
        for hdxm, d_calc_s in zip(self.fit_result.data_obj, d_calc):
            name = hdxm.name
            indices = range(hdxm.Np)
            chunks = [indices[i:i + n] for i in range(0, len(indices), n)]

            tex = []
            for chunk in chunks:
                file_name = '{}.{}'.format(str(uuid.uuid4()), extension.strip('.'))
                file_path = self._temp_dir / file_name

                fig_func = partial(_peptide_uptake_figure, fig_factory, chunk, time[0], d_calc_s, hdxm)
                self.figure_queue.append((file_path, fig_func))

                tex_func = partial(_place_figure, file_path)
                tex.append(tex_func)

            self.tex_dict['peptide_uptake'][name] = tex

    def generate_latex(self, sort_by='graphs'):  # graphs = []  #todo allow for setting which graphs to output
        if sort_by == 'graphs':
            for graph_type, state_dict in self.tex_dict.items():
                #todo map graph type to human readable section name
                with self.doc.create(pyl.Section(graph_type)):
                    for state, tex_list in state_dict.items():
                        with self.doc.create(pyl.Subsection(state)):
                            [tex_func(doc=self.doc) for tex_func in tex_list]
        else:
            raise NotImplementedError('Sorting by protein state not implemented')

    def generate_figures(self, executor='process'):
        if isinstance(executor, futures.Executor):
            exec_klass = executor
        elif executor == 'process':
            exec_klass = futures.ProcessPoolExecutor()
        elif executor == 'local':
            exec_klass = LocalThreadExecutor()
        else:
            raise ValueError("Invalid value for 'executor'")

        total = len(self.figure_queue)
        ft = [exec_klass.submit(run, item) for item in self.figure_queue]
        with tqdm(total=total, desc='Generating figures') as pbar:
            for future in futures.as_completed(ft):
                pbar.update(1)

    def generate_pdf(self, file_path, cleanup=True, **kwargs):
        defaults = {'compiler_args': ['--xelatex']}
        defaults.update(kwargs)

        self.doc.generate_pdf(file_path, **defaults)

        if cleanup:
            #try:
            self._temp_dir.clean()
            #except:



def _place_figure(file_path, width=r'\textwidth', doc=None):
    with doc.create(pyl.Figure(position='H')) as tex_fig:
        tex_fig.add_image(str(file_path), width=pyl.NoEscape(width))


def _peptide_uptake_figure(fig_factory, indices, _t, _d, hdxm):
    fig, axes = fig_factory()
    axes_iter = iter(axes)  # isnt this alreay iterable?
    for i in indices:
        ax = next(axes_iter)
        ax.plot(_t, _d[i], color='r')
        ax.scatter(hdxm.timepoints, hdxm.d_exp.iloc[i], color='k')

        start, end = hdxm.coverage.data.iloc[i][['_start', '_end']]
        ax.format(title=f'Peptide_{i}: {start} - {end}')

    for ax in axes_iter:
        ax.axis('off')
    # todo second y axis with RFU
    axes.format(xscale='log', xlabel='Time (s)', ylabel='Corrected D-uptake', xformatter='log', ylim=(0, None))

    return fig


def run(item):
    file_path, fig_func = item
    fig = fig_func()
    fig.savefig(file_path)
    plt.close(fig)


class LocalThreadExecutor(futures.Executor):

    def submit(self, f, *args, **kwargs):
        future = futures.Future()
        future.set_result(f(*args, **kwargs))
        return future

    def shutdown(self, wait=True):
        pass