import os
import tempfile
import uuid
from concurrent import futures
from functools import partial
from importlib import import_module
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import proplot as pplt
import pylatex as pyl
from tqdm.auto import tqdm

from pyhdx.plot import FitResultPlotBase, SCATTER_KWARGS

geometry_options = {
    "lmargin": "1in",
    "rmargin": "1in",
    "tmargin": "1in",
}
# pyl.NoEscape
# geometry_options = {
#     'total': pyl.NoEscape(r'{170mm, 157mm}'),
#     'left': '20mm',
#     'top': '20mm'
# }

# assuming A4 210 mm width
PAGE_WIDTH = (
    210
    - pplt.units(geometry_options["lmargin"], dest="mm")
    - pplt.units(geometry_options["rmargin"], dest="mm")
)
# PAGE_WIDTH = 170


class BaseReport(object):
    pass


class Report(BaseReport):
    def __init__(self, hdxm_set, **kwargs):
        raise NotImplementedError()


class FitReport(FitResultPlotBase):
    """
    Create .pdf output of a fit result
    """

    def __init__(
        self, fit_result, title=None, doc=None, add_date=True, temp_dir=None, **kwargs
    ):
        super().__init__(fit_result)
        self.title = title or f"Fit report"
        self.doc = doc or self._init_doc(add_date=add_date)
        self._temp_dir = temp_dir or self.make_temp_dir()
        self._temp_dir = Path(self._temp_dir)

        self.figure_queue = []
        self.tex_dict = (
            {}
        )  # dictionary gathering lists of partial functions which when executed generate the tex output
        self._figure_number = 0  # todo automate

    def make_temp_dir(self):
        # todo pathlib
        _tmp_path = os.path.abspath(os.path.join(tempfile.gettempdir(), str(id(self))))

        if not os.path.exists(_tmp_path):
            os.makedirs(_tmp_path)
        return _tmp_path

    def _init_doc(self, add_date=True):
        doc = pyl.Document(geometry_options=geometry_options)
        doc.packages.append(pyl.Package("float"))
        doc.packages.append(pyl.Package("hyperref"))

        doc.preamble.append(pyl.Command("title", self.title))
        if add_date:
            doc.preamble.append(pyl.Command("date", pyl.NoEscape(r"\today")))
        else:
            doc.preamble.append(pyl.Command("date", pyl.NoEscape(r"")))
        doc.append(pyl.NoEscape(r"\maketitle"))
        doc.append(pyl.NewPage())
        doc.append(pyl.Command("tableofcontents"))
        doc.append(pyl.NewPage())

        return doc

    # def _save_fig(self, fig, *args, extension='pdf', **kwargs):
    #     filename = '{}.{}'.format(str(uuid.uuid4()), extension.strip('.'))
    #     filepath = os.path.join(self._temp_dir, filename)
    #     fig.savefig(filepath, *args, **kwargs)
    #     return filepath

    def reset_doc(self, add_date=True):
        self.doc = self._init_doc(add_date=add_date)

    def add_standard_figure(self, name, **kwargs):
        extension = ".pdf"
        self.tex_dict[name] = {}

        module = import_module("pyhdx.plot")
        f = getattr(module, name)
        arg_dict = self._get_arg(name)
        width = kwargs.pop("width", PAGE_WIDTH)

        for args_name, arg in arg_dict.items():
            fig_func = partial(
                f, arg, width=width, **kwargs
            )  # todo perhaps something like fig = lazy(func(args, **kwargs))?
            file_name = "{}.{}".format(str(uuid.uuid4()), extension.strip("."))
            file_path = self._temp_dir / file_name

            self.figure_queue.append((file_path, fig_func))

            tex_func = partial(_place_figure, file_path)
            self.tex_dict[name][args_name] = [tex_func]

    # def _get_args(self, plot_func_name):
    #     #Add _figure suffix if not present
    #     if not plot_func_name.endswith('_figure'):
    #         plot_func_name += '_figure'
    #
    #     if plot_func_name == 'peptide_coverage_figure':
    #         return {hdxm.name: [hdxm.data] for hdxm in self.fit_result.hdxm_set.hdxm_list}
    #     elif plot_func_name == 'residue_time_scatter_figure':
    #         return {hdxm.name: [hdxm] for hdxm in self.fit_result.hdxm_set.hdxm_list}
    #     elif plot_func_name == 'residue_scatter_figure':
    #         return {'All states': [self.fit_result.hdxm_set]}
    #     elif plot_func_name == 'dG_scatter_figure':
    #         return {'All states': [self.fit_result.output]}
    #     elif plot_func_name == 'ddG_scatter_figure':
    #         return {'All states': [self.fit_result.output.df]}  # Todo change protein object to dataframe!
    #     elif plot_func_name == 'linear_bars_figure':
    #         return {'All states': [self.fit_result.output.df]}
    #     elif plot_func_name == 'rainbowclouds_figure':
    #         return {'All states': [self.fit_result.output.df]}
    #     else:
    #         raise ValueError(f"Unknown plot function {plot_func_name!r}")

    def add_peptide_uptake_curves(self, layout=(5, 4), time_axis=None):
        extension = ".pdf"
        self.tex_dict["peptide_uptake"] = {}

        nrows, ncols = layout
        n = nrows * ncols
        time = time_axis or self.get_fit_timepoints()
        if time.ndim == 1:
            time = np.tile(
                time, (len(self.fit_result), 1)
            )  # todo move shape change to FitResult object

        d_calc = self.fit_result(time)  # Ns x Np x Nt

        fig_factory = partial(
            pplt.subplots,
            ncols=ncols,
            nrows=nrows,
            sharex=1,
            sharey=1,
            width=f"{PAGE_WIDTH}mm",
        )

        # iterate over samples
        for hdxm, d_calc_s in zip(self.fit_result.hdxm_set, d_calc):
            name = hdxm.name
            indices = range(hdxm.Np)
            chunks = [indices[i : i + n] for i in range(0, len(indices), n)]

            tex = []
            for chunk in chunks:
                file_name = "{}.{}".format(str(uuid.uuid4()), extension.strip("."))
                file_path = self._temp_dir / file_name

                fig_func = partial(
                    _peptide_uptake_figure, fig_factory, chunk, time[0], d_calc_s, hdxm
                )
                self.figure_queue.append((file_path, fig_func))

                tex_func = partial(_place_figure, file_path)
                tex.append(tex_func)

            self.tex_dict["peptide_uptake"][name] = tex

    def generate_latex(
        self, sort_by="graphs"
    ):  # graphs = []  #todo allow for setting which graphs to output
        if sort_by == "graphs":
            for graph_type, state_dict in self.tex_dict.items():
                # todo map graph type to human readable section name
                with self.doc.create(pyl.Section(graph_type)):
                    for state, tex_list in state_dict.items():
                        with self.doc.create(pyl.Subsection(state)):
                            [tex_func(doc=self.doc) for tex_func in tex_list]
        else:
            raise NotImplementedError("Sorting by protein state not implemented")

    def generate_figures(self, executor="local"):
        if isinstance(executor, futures.Executor):
            exec_klass = executor
        elif executor == "process":
            exec_klass = futures.ProcessPoolExecutor()
        elif executor == "local":
            exec_klass = LocalThreadExecutor()
        else:
            raise ValueError("Invalid value for 'executor'")

        total = len(self.figure_queue)
        ft = [exec_klass.submit(run, item) for item in self.figure_queue]
        with tqdm(total=total, desc="Generating figures") as pbar:
            for future in futures.as_completed(ft):
                pbar.update(1)

    def generate_pdf(self, file_path, cleanup=True, **kwargs):
        defaults = {"compiler_args": ["--xelatex"]}
        defaults.update(kwargs)

        self.doc.generate_pdf(file_path, **defaults)
        #
        # if cleanup:
        #     #try:
        #     self._temp_dir.clean()
        #     #except:


def _place_figure(file_path, width=r"\textwidth", doc=None):
    with doc.create(pyl.Figure(position="H")) as tex_fig:
        tex_fig.add_image(str(file_path), width=pyl.NoEscape(width))


def _peptide_uptake_figure(fig_factory, indices, _t, _d, hdxm):
    fig, axes = fig_factory()
    axes_iter = iter(axes)
    for i in indices:
        ax = next(axes_iter)
        ax.plot(_t, _d[i], color="r")
        ax.scatter(hdxm.timepoints, hdxm.d_exp.iloc[i], color="k", **SCATTER_KWARGS)

        start, end = hdxm.coverage.data.iloc[i][["_start", "_end"]]
        ax.format(title=f"Peptide_{i}: {start} - {end}")

    for ax in axes_iter:
        ax.axis("off")
    # todo second y axis with RFU
    axes.format(
        xscale="log",
        xlabel="Time (s)",
        ylabel="Corrected D-uptake",
        xformatter="log",
        ylim=(0, None),
    )

    return fig


def run(item):
    file_path, fig_func = item
    fig = fig_func()
    if not isinstance(fig, plt.Figure):
        fig = fig[0]
    fig.savefig(file_path)
    plt.close(fig)


class LocalThreadExecutor(futures.Executor):
    def submit(self, f, *args, **kwargs):
        future = futures.Future()
        future.set_result(f(*args, **kwargs))
        return future

    def shutdown(self, wait=True):
        pass
