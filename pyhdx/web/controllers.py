import asyncio
import sys
import uuid
import warnings
import zipfile
from datetime import datetime
from functools import partial
from io import StringIO, BytesIO

import colorcet
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import panel as pn
import param
import yaml
from distributed import Client
from panel.io.server import async_execute
from proplot import to_hex
from skimage.filters import threshold_multiotsu

from pyhdx.batch_processing import StateParser
from pyhdx.config import cfg
from pyhdx.fileIO import read_dynamx, csv_to_dataframe, dataframe_to_stringio
from pyhdx.fitting import (
    fit_rates_weighted_average,
    fit_rates_half_time_interpolate,
    get_bounds,
    fit_gibbs_global,
    fit_gibbs_global_batch,
    PATIENCE,
    STOP_LOSS,
    EPOCHS,
    R1,
    R2,
    optimizer_defaults,
    RatesFitResult,
)
from pyhdx.models import PeptideMasterTable, HDXMeasurement, array_intersection
from pyhdx.fitting_torch import TorchFitResultSet
from pyhdx.plot import (
    dG_scatter_figure,
    ddG_scatter_figure,
    linear_bars_figure,
    rainbowclouds_figure,
    CMAP_NORM_DEFAULTS,
)
from pyhdx.support import (
    series_to_pymol,
    apply_cmap,
    multiindex_astype,
    multiindex_set_categories,
)
from pyhdx.web.base import ControlPanel, DEFAULT_CLASS_COLORS
from pyhdx.web.opts import CmapOpts
from pyhdx.web.transforms import CrossSectionTransform
from pyhdx.web.utils import fix_multiindex_dtypes
from pyhdx.web.widgets import ASyncProgressBar


def blocking_function(duration):
    import time

    time.sleep(duration)

    df = pd.DataFrame(
        {
            "x": np.random.normal(loc=3, scale=2, size=100),
            "y": np.random.normal(loc=2, scale=0.3, size=100),
        }
    )

    return df


class AsyncControlPanel(ControlPanel):
    _type = "async"

    btn = param.Action(lambda self: self.do_stuff())

    print = param.Action(lambda self: self.print_stuff())

    # async_do = param.Action(lambda self:)

    value = param.String()

    slider = param.Number(2.4, bounds=(1.0, 4.0))

    field = param.String(
        default="init_value", doc="what is the value of this field during async?"
    )

    @property
    def src(self):
        return self.sources["main"]

    def make_dict(self):
        widgets = self.generate_widgets()
        widgets["pbar"] = ASyncProgressBar()

        return widgets

    async def work_func(self):
        print("start")
        print(self.field)
        name = self.field  # is indeed stored locally
        scheduler_address = cfg.get("cluster", "scheduler_address")
        async with Client(scheduler_address, asynchronous=True) as client:
            futures = []
            for i in range(10):
                duration = (i + np.random.rand()) / 3.0
                future = client.submit(blocking_function, duration)
                futures.append(future)

            await self.widgets["pbar"].run(futures)

            results = await asyncio.gather(*futures)

        result = pd.concat(results)

        print(result)
        print(self.field)
        print(name)

        self.src.add_table("test_data", result)
        self.src.updated = True

        # self.src.param.trigger('updated')

    def do_stuff(self):

        async_execute(self.work_func)

    @param.depends("slider", watch=True)
    def _slider_updated(self):
        self.value = str(self.slider)

    def print_stuff(self):
        print(self.parent.executor)
        df = self.src.tables["test_data"]
        print()


class DevTestControl(ControlPanel):

    header = "Debug"

    _type = "dev"

    debug_btn = param.Action(lambda self: self._action_debug(), label="Debug")

    test_btn = param.Action(lambda self: self._action_test(), label="Test")

    def __init__(self, parent, **params):
        super().__init__(parent, **params)

    @property
    def src(self):
        return self.sources["main"]

    def _action_debug(self):
        transforms = self.transforms
        sources = self.sources
        views = self.views
        opts = self.opts

        tables = self.sources["main"].tables
        rfus = tables["rfu_residues"]
        print(rfus)
        print(rfus.columns.dtypes)

        c = self.parent.control_panels

        drfu = tables.get("drfu_comparison")
        if drfu is not None:
            print(drfu.columns.dtypes)
            print(drfu)

        drfu_selected = self.transforms["drfu_comparison_select"].get()
        print("break")

    def _action_test(self):
        pdbe_view = self.views["protein"]
        pdbe_view.pdbe.test = not pdbe_view.pdbe.test

    @property
    def _layout(self):
        return [
            ("self", None),
            # ('transforms.protein_select', None),
            # ('transforms.table_2_select', None)
        ]


class MWEControl(ControlPanel):
    """Temporary controller for use in the MWE app"""

    _type = "mwe"

    new_data = param.Action(lambda self: self._action_new_data())

    @property
    def src(self):
        return self.sources["main"]

    def _action_new_data(self):
        df = pd.DataFrame(
            {
                "x": np.random.normal(loc=3, scale=2, size=100),
                "y": np.random.normal(loc=2, scale=0.3, size=100),
                "yerr": np.random.uniform(0.2, 1.0, size=100),
            }
        )

        self.src.tables["test_data"] = df
        self.src.param.trigger("updated")


class PyHDXControlPanel(ControlPanel):
    @property
    def src(self):
        return self.parent.sources["main"]


class PeptideFileInputControl(PyHDXControlPanel):
    """
    This controller allows users to input .csv file (Currently only DynamX format) of 'state' peptide uptake data.
    Users can then choose how to correct for back-exchange and which 'state' and exposure times should be used for
    analysis.

    """

    _type = "peptide_file_input"

    header = "Peptide Input"

    input_mode = param.Selector(default="Manual", objects=["Manual", "Batch"])

    input_files_label = param.String("Input files:")

    input_files = param.List(
        doc="HDX input files. Currently only support DynamX format"
    )

    batch_file_label = param.String("Batch file (yaml)")

    batch_file = param.Parameter(doc="Batch file input:")

    be_mode = param.Selector(
        doc="Select method of back exchange correction",
        label="Back exchange correction method",
        objects=["FD Sample", "Flat percentage"],
    )

    fd_state = param.Selector(doc="State used to normalize uptake", label="FD State")

    fd_exposure = param.Selector(
        doc="Exposure used to normalize uptake", label="FD Exposure"
    )

    be_percent = param.Number(
        28.0,
        bounds=(0, 100),
        doc="Global percentage of back-exchange",
        label="Back exchange percentage",
    )

    exp_state = param.Selector(
        doc="State for selected experiment", label="Experiment State"
    )

    exp_exposures = param.ListSelector(
        default=[],
        objects=[""],
        label="Experiment Exposures",
        doc="Selected exposure time to use",
    )

    drop_first = param.Integer(
        1, bounds=(0, None), doc="Select the number of N-terminal residues to ignore."
    )

    d_percentage = param.Number(
        90.0,
        bounds=(0, 100),
        doc="Percentage of deuterium in the labelling buffer",
        label="Deuterium percentage",
    )

    temperature = param.Number(
        293.15,
        bounds=(0, 373.15),
        doc="Temperature of the D-labelling reaction",
        label="Temperature (K)",
    )

    pH = param.Number(
        7.5,
        doc="pH of the D-labelling reaction, as read from pH meter",
        label="pH read",
    )

    n_term = param.Integer(
        1,
        doc="Index of the n terminal residue in the protein. Can be set to negative values to "
        "accommodate for purification tags. Used in the determination of intrinsic rate of exchange",
    )

    c_term = param.Integer(
        0,
        bounds=(0, None),
        doc="Index of the c terminal residue in the protein. Used for generating pymol export script"
        "and determination of intrinsic rate of exchange for the C-terminal residue",
    )

    sequence = param.String("", doc="Optional FASTA protein sequence")

    dataset_name = param.String()

    add_dataset_button = param.Action(
        lambda self: self._action_add_dataset(),
        label="Add measurement(s)",
        doc="Parse selected peptides for further analysis and apply back-exchange correction",
    )

    hdxm_list = param.ListSelector(
        label="HDX Measurements", doc="Lists added HDX-MS measurements", constant=True
    )

    def __init__(self, parent, **params):
        super(PeptideFileInputControl, self).__init__(
            parent, _excluded=["be_percent", "batch_file", "batch_file_label"], **params
        )
        self.src.param.watch(self._hdxm_objects_updated, ["hdxm_objects"])
        self.update_box()

        self._df = None  # Dataframe with raw input data

    def make_dict(self):
        text_area = pn.widgets.TextAreaInput(
            name="Sequence (optional)",
            placeholder="Enter sequence in FASTA format",
            max_length=10000,
            width=300,
            height=100,
            height_policy="fixed",
            width_policy="fixed",
        )
        return self.generate_widgets(
            input_files_label=pn.widgets.StaticText(value=self.input_files_label),
            input_files=pn.widgets.FileInput(multiple=True, name="Input files"),
            batch_file_label=pn.widgets.StaticText(value=self.batch_file_label),
            batch_file=pn.widgets.FileInput(name="Batch yaml file", accept=".yaml"),
            temperature=pn.widgets.FloatInput,
            be_percent=pn.widgets.FloatInput,
            d_percentage=pn.widgets.FloatInput,
            sequence=text_area,
        )

    def make_list(self):
        excluded = ["be_percent"]
        widget_list = [
            widget for name, widget, in self.widget_dict.items() if name not in excluded
        ]

        return widget_list

    @param.depends("be_mode", "input_mode", watch=True)
    def _update_mode(self):
        excluded = set()
        if self.be_mode == "FD Sample":
            excluded |= {"be_percent"}
        elif self.be_mode == "Flat percentage":
            excluded |= {"fd_state", "fd_exposure"}

        if self.input_mode == "Manual":
            excluded |= {"batch_file", "batch_file_label"}
        elif self.input_mode == "Batch":
            excluded |= {
                "be_mode",
                "fd_state",
                "fd_exposure",
                "be_percent",
                "exp_state",
                "exp_exposures",
                "drop_first",
                "d_percentage",
                "temperature",
                "pH",
                "n_term",
                "c_term",
                "sequence",
                "dataset_name",
            }

        self._excluded = list(excluded)

        self.update_box()

    @param.depends("input_files", watch=True)
    def _read_files(self):
        if self.input_files:
            combined_df = read_dynamx(
                *[
                    StringIO(byte_content.decode("UTF-8"))
                    for byte_content in self.input_files
                ]
            )
            self._df = combined_df

            self.parent.logger.info(
                f'Loaded {len(self.input_files)} file{"s" if len(self.input_files) > 1 else ""} with a total '
                f"of {len(self._df)} peptides"
            )

        else:
            self._df = None

        self._update_fd_state()
        self._update_fd_exposure()
        self._update_exp_state()
        self._update_exp_exposure()

    def _update_fd_state(self):
        if self._df is not None:
            states = list(self._df["state"].unique())
            self.param["fd_state"].objects = states
            self.fd_state = states[0]
        else:
            self.param["fd_state"].objects = []

    @param.depends("fd_state", watch=True)
    def _update_fd_exposure(self):
        if self._df is not None:
            fd_entries = self._df[self._df["state"] == self.fd_state]
            exposures = list(np.unique(fd_entries["exposure"]))
        else:
            exposures = []
        self.param["fd_exposure"].objects = exposures
        if exposures:
            self.fd_exposure = exposures[0]

    @param.depends("fd_state", "fd_exposure", watch=True)
    def _update_exp_state(self):
        if self._df is not None:
            # Booleans of data entries which are in the selected control
            control_bools = np.logical_and(
                self._df["state"] == self.fd_state,
                self._df["exposure"] == self.fd_exposure,
            )

            control_data = self._df[control_bools].to_records()
            other_data = self._df[~control_bools].to_records()

            intersection = array_intersection(
                [control_data, other_data], fields=["start", "end"]
            )  # sequence?
            states = list(np.unique(intersection[1]["state"]))
        else:
            states = []

        self.param["exp_state"].objects = states
        if states:
            self.exp_state = states[0] if not self.exp_state else self.exp_state

    @param.depends("exp_state", watch=True)
    def _update_exp_exposure(self):
        if self._df is not None:
            exp_entries = self._df[self._df["state"] == self.exp_state]
            exposures = list(np.unique(exp_entries["exposure"]))
            exposures.sort()
        else:
            exposures = []

        self.param["exp_exposures"].objects = exposures
        self.exp_exposures = [e for e in exposures if e != 0.0]

        if (
            not self.dataset_name
            or self.dataset_name in self.param["exp_state"].objects
        ):
            self.dataset_name = self.exp_state

        if not self.c_term and exposures:
            self.c_term = int(np.max(exp_entries["end"]))

    def _hdxm_objects_updated(self, events):
        # Update datasets widget as datasets on parents change
        objects = list(self.src.hdxm_objects.keys())
        self.param["hdxm_list"].objects = objects

    def _action_add_dataset(self):
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.HDXMeasurement`"""
        if self.input_mode == "Manual":
            self._add_dataset_manual()
        elif self.input_mode == "Batch":
            self._add_dataset_batch()

    def _add_dataset_batch(self):
        if self._df is None:
            self.parent.logger.info("No data loaded")
            return
        if self.src.hdxm_objects:
            self.parent.logger.info(
                "Can only batch load data when no data was previously loaded"
            )
            return

        yaml_dict = yaml.safe_load(self.batch_file.decode("UTF-8"))
        ios = {
            name: StringIO(byte_content.decode("UTF-8"))
            for name, byte_content in zip(
                self.widgets["input_files"].filename, self.input_files
            )
        }
        filters = [lambda df: df.query("exposure > 0")]

        parser = StateParser(yaml_dict, data_src=ios, data_filters=filters)

        for state in yaml_dict.keys():
            hdxm = parser.load_hdxm(state, name=state)
            self.src.add(hdxm, state)
            self.parent.logger.info(
                f"Loaded dataset {state} with experiment state {hdxm.state} "
                f"({len(hdxm)} timepoints, {len(hdxm.coverage)} peptides each)"
            )
            self.parent.logger.info(
                f"Average coverage: {hdxm.coverage.percent_coverage:.3}%, "
                f"Redundancy: {hdxm.coverage.redundancy:.2}"
            )

    def _add_dataset_manual(self):

        if self._df is None:
            self.parent.logger.info("No data loaded")
            return
        elif self.dataset_name in self.src.hdxm_objects.keys():
            self.parent.logger.info(f"Dataset name {self.dataset_name} already in use")
            return

        peptides = PeptideMasterTable(
            self._df, d_percentage=self.d_percentage, drop_first=self.drop_first,
        )
        if self.be_mode == "FD Sample":
            control_0 = None  # = (self.zero_state, self.zero_exposure) if self.zero_state != 'None' else None
            peptides.set_control((self.fd_state, self.fd_exposure), control_0=control_0)
        elif self.be_mode == "Flat percentage":
            # todo @tejas: Add test
            peptides.set_backexchange(self.be_percent)

        data = peptides.get_state(self.exp_state)
        exp_bools = data["exposure"].isin(self.exp_exposures)
        data = data[exp_bools]

        # todo temperature ph kwarg for series
        hdxm = HDXMeasurement(
            data,
            c_term=self.c_term,
            n_term=self.n_term,
            sequence=self.sequence,
            name=self.dataset_name,
            temperature=self.temperature,
            pH=self.pH,
        )

        self.src.add(hdxm, self.dataset_name)
        self.parent.logger.info(
            f"Loaded dataset {self.dataset_name} with experiment state {self.exp_state} "
            f"({len(hdxm)} timepoints, {len(hdxm.coverage)} peptides each)"
        )
        self.parent.logger.info(
            f"Average coverage: {hdxm.coverage.percent_coverage:.3}%, "
            f"Redundancy: {hdxm.coverage.redundancy:.2}"
        )

    def _action_remove_datasets(self):
        raise NotImplementedError("Removing datasets not implemented")
        for name in self.hdxm_list:
            self.parent.datasets.pop(name)

        self.parent.param.trigger(
            "datasets"
        )  # Manual trigger as key assignment does not trigger the param


class PeptideRFUFileInputControl(PyHDXControlPanel):
    """
    This controller allows users to input .csv file (Currently only DynamX format) of 'state' peptide uptake data.
    Users can then choose how to correct for back-exchange and which 'state' and exposure times should be used for
    analysis.

    """

    _type = "peptide_rfu_file_input"

    header = "Peptide Input"

    input_files = param.List()

    fd_state = param.Selector(doc="State used to normalize uptake", label="FD State")

    fd_exposure = param.Selector(
        doc="Exposure used to normalize uptake", label="FD Exposure"
    )

    nd_state = param.Selector(doc="State used to normalize uptake", label="ND State")

    nd_exposure = param.Selector(
        doc="Exposure used to normalize uptake", label="ND Exposure"
    )

    exp_state = param.Selector(
        doc="State for selected experiment", label="Experiment State"
    )

    exp_exposures = param.ListSelector(
        default=[],
        objects=[""],
        label="Experiment Exposures",
        doc="Selected exposure time to use",
    )

    d_percentage = param.Number(
        95.0,
        bounds=(0, 100),
        doc="Percentage of deuterium in the labelling buffer",
        label="Deuterium percentage",
    )

    drop_first = param.Integer(
        1, bounds=(0, None), doc="Select the number of N-terminal residues to ignore."
    )

    n_term = param.Integer(
        1,
        doc="Index of the n terminal residue in the protein. Can be set to negative values to "
        "accommodate for purification tags. Used in the determination of intrinsic rate of exchange",
    )

    c_term = param.Integer(
        0,
        bounds=(0, None),
        doc="Index of the c terminal residue in the protein. Used for generating pymol export script"
        "and determination of intrinsic rate of exchange for the C-terminal residue",
    )

    sequence = param.String("", doc="Optional FASTA protein sequence")

    dataset_name = param.String()

    add_dataset_button = param.Action(
        lambda self: self._action_add_dataset(),
        label="Add measurement",
        doc="Parse selected peptides for further analysis and apply back-exchange correction",
    )

    hdxm_list = param.ListSelector(
        label="HDX Measurements", doc="Lists added HDX-MS measurements", constant=True
    )

    def __init__(self, parent, **params):
        super(PeptideRFUFileInputControl, self).__init__(
            parent, _excluded=["be_percent"], **params
        )
        self.src.param.watch(self._hdxm_objects_updated, ["hdxm_objects"])
        self.update_box()

        self._df = None  # Numpy array with raw input data (or is pd.Dataframe?)

    def make_dict(self):
        text_area = pn.widgets.TextAreaInput(
            name="Sequence (optional)",
            placeholder="Enter sequence in FASTA format",
            max_length=10000,
            width=300,
            height=100,
            height_policy="fixed",
            width_policy="fixed",
        )
        return self.generate_widgets(
            input_files=pn.widgets.FileInput(multiple=True, name="Input files"),
            temperature=pn.widgets.FloatInput,
            d_percentage=pn.widgets.FloatInput,
            # fd_percentage=pn.widgets.FloatInput,
            sequence=text_area,
        )

    @param.depends("input_files", watch=True)
    def _read_files(self):
        if self.input_files:
            combined_df = read_dynamx(
                *[
                    StringIO(byte_content.decode("UTF-8"))
                    for byte_content in self.input_files
                ]
            )
            self._df = combined_df

            self.parent.logger.info(
                f'Loaded {len(self.input_files)} file{"s" if len(self.input_files) > 1 else ""} with a total '
                f"of {len(self._df)} peptides"
            )
        else:
            self._df = None

        self._update_fd_state()
        self._update_fd_exposure()
        self._update_nd_state()
        self._update_nd_exposure()
        self._update_exp_state_fd()
        self._update_exp_state_nd()
        self._update_exp_exposure()

    def _update_nd_state(self):
        if self._df is not None:
            states = list(self._df["state"].unique())
            self.param["nd_state"].objects = states
            self.nd_state = states[0]
        else:
            self.param["fd_state"].objects = []

    @param.depends("nd_state", watch=True)
    def _update_nd_exposure(self):
        if self._df is not None:
            nd_entries = self._df[self._df["state"] == self.nd_state]
            exposures = list(np.unique(nd_entries["exposure"]))
        else:
            exposures = []
        self.param["nd_exposure"].objects = exposures
        if exposures:
            self.nd_exposure = exposures[0]

    def _update_fd_state(self):
        if self._df is not None:
            states = list(self._df["state"].unique())
            self.param["fd_state"].objects = states
            self.fd_state = states[0]
        else:
            self.param["fd_state"].objects = []

    @param.depends("fd_state", watch=True)
    def _update_fd_exposure(self):
        if self._df is not None:
            fd_entries = self._df[self._df["state"] == self.fd_state]
            exposures = list(np.unique(fd_entries["exposure"]))
        else:
            exposures = []
        self.param["fd_exposure"].objects = exposures
        if exposures:
            self.fd_exposure = exposures[0]

    @param.depends("fd_state", "fd_exposure", watch=True)
    def _update_exp_state_fd(self):
        if self._df is not None:
            # Booleans of data entries which are in the selected control
            fd_bools = np.logical_and(
                self._df["state"] == self.fd_state,
                self._df["exposure"] == self.fd_exposure,
            )

            control_data = self._df[fd_bools].to_records()
            other_data = self._df[~fd_bools].to_records()

            intersection = array_intersection(
                [control_data, other_data], fields=["start", "end"]
            )  # sequence?
            states = list(np.unique(intersection[1]["state"]))
        else:
            states = []

        self.param["exp_state"].objects = states
        if states:
            self.exp_state = states[0] if not self.exp_state else self.exp_state

    @param.depends("nd_state", "nd_exposure", watch=True)
    def _update_exp_state_nd(self):
        if self._df is not None:
            # Booleans of data entries which are in the selected control
            fd_bools = np.logical_and(
                self._df["state"] == self.nd_state,
                self._df["exposure"] == self.nd_exposure,
            )

            control_data = self._df[fd_bools].to_records()
            other_data = self._df[~fd_bools].to_records()

            intersection = array_intersection(
                [control_data, other_data], fields=["start", "end"]
            )  # sequence?
            states = list(np.unique(intersection[1]["state"]))
        else:
            states = []

        self.param["exp_state"].objects = states
        if states:
            self.exp_state = states[0] if not self.exp_state else self.exp_state

    @param.depends("exp_state", watch=True)
    def _update_exp_exposure(self):
        if self._df is not None:
            exp_entries = self._df[self._df["state"] == self.exp_state]
            exposures = list(np.unique(exp_entries["exposure"]))
            exposures.sort()
        else:
            exposures = []

        self.param["exp_exposures"].objects = exposures
        self.exp_exposures = [e for e in exposures if e != 0.0]

        if (
            not self.dataset_name
            or self.dataset_name in self.param["exp_state"].objects
        ):
            self.dataset_name = self.exp_state

        if not self.c_term and exposures:
            self.c_term = int(np.max(exp_entries["end"]))

    def _hdxm_objects_updated(self, events):
        # Update datasets widget as datasets on parents change
        objects = list(self.src.hdxm_objects.keys())
        self.param["hdxm_list"].objects = objects

    def _action_add_dataset(self):
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.HDXMeasurement`"""

        if self._df is None:
            self.parent.logger.info("No data loaded")
            return
        elif self.dataset_name in self.src.hdxm_objects.keys():
            self.parent.logger.info(f"Dataset name {self.dataset_name} already in use")
            return

        peptides = PeptideMasterTable(
            self._df, d_percentage=self.d_percentage, drop_first=self.drop_first,
        )

        peptides.set_control(
            (self.fd_state, self.fd_exposure),
            control_0=(self.nd_state, self.nd_exposure),
        )

        data = peptides.get_state(self.exp_state)
        exp_bools = data["exposure"].isin(self.exp_exposures)
        data = data[exp_bools]

        # todo temperature ph kwarg for series
        hdxm = HDXMeasurement(
            data,
            c_term=self.c_term,
            n_term=self.n_term,
            sequence=self.sequence,
            name=self.dataset_name,
        )

        self.src.add(hdxm, self.dataset_name)
        self.parent.logger.info(
            f"Loaded dataset {self.dataset_name} with experiment state {self.exp_state} "
            f"({len(hdxm)} timepoints, {len(hdxm.coverage)} peptides each)"
        )
        self.parent.logger.info(
            f"Average coverage: {hdxm.coverage.percent_coverage:.3}%, "
            f"Redundancy: {hdxm.coverage.redundancy:.2}"
        )

    def _action_remove_datasets(self):
        raise NotImplementedError("Removing datasets not implemented")
        for name in self.hdxm_list:
            self.parent.datasets.pop(name)

        self.parent.param.trigger(
            "datasets"
        )  # Manual trigger as key assignment does not trigger the param


class InitialGuessControl(PyHDXControlPanel):
    """
    This controller allows users to derive initial guesses for D-exchange rate from peptide uptake data.
    """

    _type = "initial_guess"

    header = "Initial Guesses"
    fitting_model = param.Selector(
        default="Half-life (λ)",
        objects=["Half-life (λ)", "Association"],
        doc="Choose method for determining initial guesses.",
    )
    dataset = param.Selector(
        default="", doc="Dataset to apply bounds to", label="Dataset (for bounds)"
    )
    global_bounds = param.Boolean(
        default=False, doc="Set bounds globally across all datasets"
    )
    lower_bound = param.Number(0.0, doc="Lower bound for association model fitting")

    upper_bound = param.Number(0.0, doc="Upper bound for association model fitting")

    guess_name = param.String(default="Guess_1", doc="Name for the initial guesses")

    _guess_names = param.List(
        [], doc="List of current and future guess names", precedence=-1
    )

    do_fit1 = param.Action(
        lambda self: self._action_fit(),
        label="Calculate Guesses",
        doc="Start initial guess fitting",
        constant=True,
    )

    bounds = param.Dict(
        {}, doc="Dictionary which stores rate fitting bounds", precedence=-1
    )

    def __init__(self, parent, **params):
        _excluded = ["lower_bound", "upper_bound", "global_bounds", "dataset"]
        super(InitialGuessControl, self).__init__(parent, _excluded=_excluded, **params)
        self.src.param.watch(
            self._parent_hdxm_objects_updated, ["hdxm_objects"]
        )  # todo refactor ( to what?)

        self.update_box()

    def make_dict(self):
        widgets = self.generate_widgets(
            lower_bound=pn.widgets.FloatInput, upper_bound=pn.widgets.FloatInput
        )

        widgets["pbar"] = ASyncProgressBar()

        return widgets

    @param.depends("fitting_model", watch=True)
    def _fitting_model_updated(self):
        if self.fitting_model == "Half-life (λ)":
            self._excluded = ["dataset", "lower_bound", "upper_bound", "global_bounds"]

        elif self.fitting_model in ["Association", "Dissociation"]:
            self._excluded = []

        self.update_box()

    @param.depends("global_bounds", watch=True)
    def _global_bounds_updated(self):
        if self.global_bounds:
            self.param["dataset"].constant = True
        else:
            self.param["dataset"].constant = False

    @param.depends("dataset", watch=True)
    def _dataset_updated(self):
        lower, upper = self.bounds[self.dataset]
        self.lower_bound = lower
        self.upper_bound = upper

    @param.depends("lower_bound", "upper_bound", watch=True)
    def _bounds_updated(self):
        if not self.global_bounds:
            self.bounds[self.dataset] = (self.lower_bound, self.upper_bound)

    def _parent_hdxm_objects_updated(self, *events):
        if len(self.src.hdxm_objects) > 0:
            self.param["do_fit1"].constant = False

        # keys to remove:
        for k in self.bounds.keys() - self.src.hdxm_objects.keys():
            self.bounds.pop(k)
        # keys to add:
        for k in self.src.hdxm_objects.keys() - self.bounds.keys():
            self.bounds[k] = get_bounds(self.src.hdxm_objects[k].timepoints)

        options = list(self.src.hdxm_objects.keys())
        self.param["dataset"].objects = options
        if not self.dataset:
            self.dataset = options[0]

    def add_fit_result(self, future):
        name = self._guess_names.pop(future.key)

        results = future.result()
        result_obj = RatesFitResult(results)
        self.src.add(result_obj, name)

        self.param["do_fit1"].constant = False
        self.widgets["do_fit1"].loading = False

    def _action_fit(self):
        # Checking if data is available, should always be the case as button is locked before
        if len(self.src.hdxm_objects) == 0:
            self.parent.logger.info("No datasets loaded")
            return

        if self.guess_name in self._guess_names:
            self.parent.logger.info(f"Guess with name {self.guess_name} already in use")
            return

        self._guess_names.append(self.guess_name)
        self.parent.logger.info("Started initial guess fit")
        self.param["do_fit1"].constant = True
        self.widgets["do_fit1"].loading = True

        num_samples = len(self.src.hdxm_objects)
        if self.fitting_model.lower() in ["association", "dissociation"]:
            loop = asyncio.get_running_loop()
            loop.create_task(self._fit_rates(self.guess_name))

        # this is practically instantaneous and does not require dask
        elif self.fitting_model == "Half-life (λ)":
            results = map(
                fit_rates_half_time_interpolate, self.src.hdxm_objects.values()
            )

            result_obj = RatesFitResult(list(results))
            self.src.add(result_obj, self.guess_name)

            self.param["do_fit1"].constant = False
            self.widgets["do_fit1"].loading = False
            self.parent.logger.info(f"Finished initial guess fit {self.guess_name}")

    async def _fit_rates(self, name):

        num_samples = len(self.src.hdxm_objects)

        scheduler_address = cfg.get("cluster", "scheduler_address")
        self.widgets["pbar"].num_tasks = num_samples
        async with Client(scheduler_address, asynchronous=True) as client:
            if self.global_bounds:
                bounds = [(self.lower_bound, self.upper_bound)] * num_samples
            else:
                bounds = self.bounds.values()
            futures = []
            for hdxm, bound in zip(self.src.hdxm_objects.values(), bounds):
                future = client.submit(
                    fit_rates_weighted_average, hdxm, bound, client="worker_client"
                )
                futures.append(future)

            await self.widgets["pbar"].run(futures)

            results = await asyncio.gather(*futures)

        result_obj = RatesFitResult(list(results))
        self.src.add(result_obj, name)

        self.param["do_fit1"].constant = False
        self.widgets["do_fit1"].loading = False
        self.parent.logger.info(f"Finished initial guess fit {name}")


class FitControl(PyHDXControlPanel):
    """
    This controller allows users to execute PyTorch fitting of the global data set.

    Currently, repeated fitting overrides the old result.
    """

    _type = "fit"

    header = "ΔG Fit"

    initial_guess = param.Selector(doc="Name of dataset to use for initial guesses.")

    guess_mode = param.Selector(
        default="One-to-one",
        objects=["One-to-one", "One-to-many"],
        doc="Use initial guesses for each protein state (one-to-one) or use one initial"
        "guess for all protein states (one-to-many)",
    )

    guess_state = param.Selector(
        doc="Which protein state to use for initial guess when using one-to-many guesses"
    )

    fit_mode = param.Selector(default="Single", objects=["Batch", "Single"])

    stop_loss = param.Number(
        STOP_LOSS,
        bounds=(0, None),
        doc="Threshold loss difference below which to stop fitting.",
    )

    stop_patience = param.Integer(
        PATIENCE,
        bounds=(1, None),
        doc="Number of epochs where stop loss should be satisfied before stopping.",
    )

    learning_rate = param.Number(
        optimizer_defaults["SGD"]["lr"],
        bounds=(0, None),
        doc="Learning rate parameter for optimization.",
    )

    momentum = param.Number(
        optimizer_defaults["SGD"]["momentum"],
        bounds=(0, None),
        doc="Stochastic Gradient Descent momentum",
    )

    nesterov = param.Boolean(
        optimizer_defaults["SGD"]["nesterov"],
        doc="Use Nesterov type of momentum for SGD",
    )

    epochs = param.Integer(
        EPOCHS, bounds=(1, None), doc="Maximum number of epochs (iterations."
    )

    r1 = param.Number(
        R1,
        bounds=(0, None),
        label="Regularizer 1 (peptide axis)",
        doc="Value of the regularizer along residue axis.",
    )

    r2 = param.Number(
        R2,
        bounds=(0, None),
        label="Regularizer 2 (sample axis)",
        doc="Value of the regularizer along sample axis.",
    )

    reference = param.Selector(
        None,
        allow_None=True,
        label="R2 reference",
        doc="Select reference state to use in batch fitting",
    )

    fit_name = param.String("Gibbs_fit_1", doc="Name for for the fit result")

    _fit_names = param.List(
        [], precedence=-1, doc="List of names of completed and running fits"
    )

    do_fit = param.Action(
        lambda self: self._action_fit(),
        constant=True,
        label="Do Fitting",
        doc="Start global fitting",
    )

    def __init__(self, parent, **params):
        self.pbar1 = ASyncProgressBar()  # tqdm?
        super(FitControl, self).__init__(parent, **params)

        self.src.param.watch(self._source_updated, ["updated"])
        self._mode_updated()  # Initialize excluded widgets
        self._current_jobs = 0
        self._max_jobs = 2  # todo config

    def make_dict(self):
        widgets = self.generate_widgets()
        # widgets['pbar_col'] = pn.layout.Column()
        widgets["pbar"] = ASyncProgressBar()
        # widgets['progress'] = CallbackProgress()

        return widgets

    def _source_updated(self, *events):
        rate_objects = list(self.src.rate_results.keys())
        if rate_objects:
            self.param["do_fit"].constant = False

        fit_objects = list(self.src.dG_fits.keys())
        self.param["initial_guess"].objects = rate_objects + fit_objects
        if not self.initial_guess and rate_objects:
            self.initial_guess = rate_objects[0]

        hdxm_objects = list(self.src.hdxm_objects.keys())
        self.param["reference"].objects = [None] + hdxm_objects
        self.param["guess_state"].objects = hdxm_objects
        if not self.guess_state and hdxm_objects:
            self.guess_state = hdxm_objects[0]

        self._mode_updated()

    @param.depends("guess_mode", "fit_mode", watch=True)
    def _mode_updated(self):
        excluded = []
        if not (self.fit_mode == "Batch" and len(self.src.hdxm_objects) > 1):
            excluded += ["r2", "reference"]
        if self.guess_mode == "One-to-one":
            excluded += ["guess_state"]
        self._excluded = excluded
        self.update_box()

    # @param.depends("fit_mode", watch=True)
    # def _fit_mode_updated(self):
    #     if self.fit_mode == "Batch" and len(self.src.hdxm_objects) > 1:
    #         # self.param['r2'].constant = False
    #         self._excluded = []
    #     else:
    #         # self.param['r2'].constant = True
    #         self._excluded = ["r2", "reference"]
    #
    #     self.update_box()

    def _action_fit(self):
        if self.fit_name in self._fit_names:
            self.parent.logger.info(
                f"Fit result with name {self.fit_name} already in use"
            )
            return
        self._fit_names.append(self.fit_name)
        self.parent.logger.info("Started PyTorch fit")

        self._current_jobs += 1
        # if self._current_jobs >= self._max_jobs:
        #     self.widgets['do_fit'].constant = True

        self.widgets["do_fit"].loading = True

        self.parent.logger.info(f"Current number of active jobs: {self._current_jobs}")

        if self.fit_mode == "Batch":
            async_execute(self._batch_fit)
        else:
            async_execute(self._single_fit)

    def get_guesses(self):
        ...

        # initial guesses are rates
        if self.initial_guess in self.src.rate_results:
            rates_df = self.src.get_table("rates")

            if self.guess_mode == "One-to-one":
                sub_df = rates_df.xs((self.initial_guess, "rate"), level=[0, 2], axis=1)
                gibbs_guess = self.src.hdx_set.guess_deltaG(sub_df)
            elif self.guess_mode == "One-to-many":
                hdxm = self.src.hdxm_objects[self.guess_state]
                rates_series = rates_df[(self.initial_guess, self.guess_state, "rate")]
                gibbs_guess = hdxm.guess_deltaG(rates_series)

        # intial guess are dG values from previous fit
        elif self.initial_guess in self.src.dG_fits:
            dG_df = self.src.get_table("dG_fits")

            if self.guess_mode == "One-to-one":
                gibbs_guess = dG_df.xs(
                    (self.initial_guess, "_dG"), level=[0, 2], axis=1
                )
            elif self.guess_mode == "One-to-many":
                gibbs_guess = dG_df[(self.initial_guess, self.guess_state, "_dG")]

        else:
            self.parent.logger.debug(f"Initial guess {self.initial_guess!r} not found")

        return gibbs_guess

    async def _single_fit(self):
        name = self.fit_name

        # data_objs = self.src.hdxm_objects.values()
        # rates_df = self.src.rate_results[self.initial_guess].output
        gibbs_guesses = (
            self.get_guesses()
        )  # returns either DataFrame or Series depending on guess mode
        futures = []

        scheduler_address = cfg.get("cluster", "scheduler_address")
        async with Client(scheduler_address, asynchronous=True) as client:
            for protein_state, hdxm in self.src.hdxm_objects.items():
                if isinstance(gibbs_guesses, pd.Series):
                    guess = gibbs_guesses
                else:
                    guess = gibbs_guesses[protein_state]

                future = client.submit(fit_gibbs_global, hdxm, guess, **self.fit_kwargs)
                futures.append(future)

            self.widgets["pbar"].num_tasks = len(futures)
            await self.widgets["pbar"].run(futures)

            results = await asyncio.gather(*futures)

        result = TorchFitResultSet(results)
        self.src.add(result, name)
        self._current_jobs -= 1
        self.widgets["pbar"].active = False
        self.widgets["do_fit"].loading = False
        self.parent.logger.info(f"Finished PyTorch fit: {name}")

    async def _batch_fit(self):
        self.widgets["pbar"].active = True
        name = self.fit_name
        hdx_set = self.src.hdx_set
        gibbs_guess = self.get_guesses()
        scheduler_address = cfg.get("cluster", "scheduler_address")
        async with Client(scheduler_address, asynchronous=True) as client:
            future = client.submit(
                fit_gibbs_global_batch, hdx_set, gibbs_guess, **self.fit_kwargs
            )
            result = await future

        self.src.add(result, name)

        self._current_jobs -= 1
        self.widgets["pbar"].active = False
        self.widgets["do_fit"].loading = False
        self.parent.logger.info(f"Finished PyTorch fit: {name}")
        self.parent.logger.info(
            f"Finished fitting in {len(result.losses)} epochs, final mean squared residuals is {result.mse_loss:.2f}"
        )
        self.parent.logger.info(
            f"Total loss: {result.total_loss:.2f}, regularization loss: {result.reg_loss:.2f} "
            f"({result.regularization_percentage:.1f}%)"
        )

    @property
    def fit_kwargs(self):
        fit_kwargs = dict(
            r1=self.r1,
            lr=self.learning_rate,
            momentum=self.momentum,
            nesterov=self.nesterov,
            epochs=self.epochs,
            patience=self.stop_patience,
            stop_loss=self.stop_loss,
        )
        # callbacks=[self.widgets['progress'].callback])
        if self.fit_mode == "Batch":
            fit_kwargs["r2"] = self.r2
            fit_kwargs["r2_reference"] = self.reference

        return fit_kwargs


class DifferentialControl(PyHDXControlPanel):
    _type = "diff"

    header = "Differential HDX"

    reference_state = param.Selector(doc="Which of the states to use as reference")

    comparison_name = param.String(
        default="comparison_1", doc="Name for the comparison table"
    )

    add_comparison = param.Action(lambda self: self._action_add_comparison())

    def __init__(self, parent, **params):
        super().__init__(parent, **params)
        self.widgets["add_comparison"].disabled = True
        self.src.param.watch(self._source_updated, "hdxm_objects")
        self._df = None
        self._source_updated()  # todo trs source does not trigger updated when init

    @property
    def _layout(self):
        layout = []
        if "ddG_fit_select" in self.transforms:
            layout.append(("transforms.ddG_fit_select", None))
        layout.append(("self", None))

        return layout

    def get(self):
        print("remove this")
        df = self.transforms["ddG_fit_select"].get()
        return df

    def _source_updated(self, *events):
        # Triggered when hdxm objects are added
        options = self.src.names
        if len(options) >= 2:
            self.widgets["add_comparison"].disabled = False

        self.param["reference_state"].objects = options
        if self.reference_state is None and options:
            self.reference_state = options[0]

    def _action_add_comparison(self):
        current_df = self.src.get_table("drfu_comparison")
        if (
            current_df is not None
            and self.comparison_name in current_df.columns.get_level_values(level=0)
        ):
            self.parent.logger.info(
                f"Comparison name {self.comparison_name!r} already exists"
            )
            return

        # RFU only app has no dGs,
        if "ddG_fit_select" in self.transforms:
            self.add_ddG_comparison()
        self.add_drfu_comparison()

        self.parent.logger.info(
            f"Successfully added comparison set {self.comparison_name!r}"
        )
        self.src.updated = True

    def add_ddG_comparison(self):
        dG_df = self.transforms["ddG_fit_select"].get()
        if dG_df is None:
            return

        reference = dG_df[self.reference_state]["dG"]
        test = dG_df.xs("dG", axis=1, level=1).drop(self.reference_state, axis=1)
        # todo repeated code in plot.ddG_scatter_figure
        ddG = test.subtract(reference, axis=0)

        names = ["comparison_name", "comparison_state", "quantity"]
        columns = pd.MultiIndex.from_product(
            [[self.comparison_name], ddG.columns, ["ddG"]], names=names
        )
        ddG.columns = columns

        cov_ref = dG_df[self.reference_state, "covariance"] ** 2
        cov_test = (
            dG_df.xs("covariance", axis=1, level=1).drop(self.reference_state, axis=1)
            ** 2
        )
        cov = cov_test.add(cov_ref, axis=0).pow(0.5)
        columns = pd.MultiIndex.from_product(
            [[self.comparison_name], cov.columns, ["covariance"]], names=names
        )
        cov.columns = columns

        combined = pd.concat([ddG, cov], axis=1)

        categories = list(combined.columns.unique(level=1))
        combined.columns = multiindex_astype(combined.columns, 1, "category")
        combined.columns = multiindex_set_categories(
            combined.columns, 1, categories, ordered=True
        )

        self.src._add_table(combined, "ddG_comparison")

        # self.parent.sources['main'].param.trigger('tables')  #todo check/remove tables trigger

    def add_drfu_comparison(self):
        # TODO adapt
        # current_df = self.src.get_table('ddG_comparison')
        # if current_df is not None and self.comparison_name in current_df.columns.get_level_values(level=0):
        #     self.parent.logger.info(f"Comparison name {self.comparison_name!r} already exists")
        #     return

        rfu_df = self.src.get_table("rfu_residues")

        reference = rfu_df[self.reference_state]
        test = (
            rfu_df.drop(self.reference_state, axis=1)
            .reorder_levels(["exposure", "state", "quantity"], axis=1)
            .sort_index(axis=1, level=0)
        )

        test = test.sort_index(axis=1, level=0)

        drfu = (
            test.sub(reference, axis="columns")
            .reorder_levels(["state", "exposure", "quantity"], axis=1)
            .sort_index(axis=1)
            .dropna(how="all", axis=1)
        )

        # Expand multiindex level
        tuples = [(self.comparison_name, *tup[:-1], "drfu") for tup in drfu.columns]
        drfu.columns = pd.MultiIndex.from_tuples(
            tuples,
            names=["comparison_name", "comparison_state", "exposure", "quantity"],
        )

        # Set the 'comparison_state' level back to categorical
        categories = list(drfu.columns.unique(level=1))
        drfu.columns = multiindex_astype(drfu.columns, 1, "category")
        drfu.columns = multiindex_set_categories(
            drfu.columns, 1, categories, ordered=True
        )

        # TODO should be public
        self.src._add_table(drfu, "drfu_comparison")


class ColorTransformControl(PyHDXControlPanel):
    """
    This controller allows users classify 'mapping' datasets and assign them colors.

    Coloring can be either in discrete categories or as a continuous custom color map.
    """

    _type = "color_transform"

    header = "Color Transform"

    # todo unify name for target field (target_data set)
    # When coupling param with the same name together there should be an option to exclude this behaviour
    quantity = param.Selector(
        label="Target Quantity"
    )  # todo refactor cmapopt / color transform??
    # fit_ID = param.Selector()  # generalize selecting widgets based on selected table
    # quantity = param.Selector(label='Quantity')  # this is the lowest-level quantity of the multiindex df (filter??)

    current_color_transform = param.String()

    mode = param.Selector(
        default="Colormap",
        objects=["Colormap", "Continuous", "Discrete"],
        doc="Choose color mode (interpolation between selected colors).",
    )
    num_colors = param.Integer(
        2,
        bounds=(1, 10),
        label="Number of colours",
        doc="Number of classification colors.",
    )
    library = param.Selector(
        default="pyhdx_default",
        objects=["pyhdx_default", "user_defined", "matplotlib", "colorcet"],
    )
    colormap = param.Selector()
    otsu_thd = param.Action(
        lambda self: self._action_otsu(),
        label="Otsu",
        doc="Automatically perform thresholding based on Otsu's method.",
    )
    linear_thd = param.Action(
        lambda self: self._action_linear(),
        label="Linear",
        doc="Automatically perform thresholding by creating equally spaced sections.",
    )
    # log_space = param.Boolean(False,
    #                          doc='Boolean to set whether to apply colors in log space or not.')
    # apply = param.Action(lambda self: self._action_apply())
    no_coverage = param.Color(
        default="#8c8c8c", doc="Color to use for regions of no coverage"
    )

    live_preview = param.Boolean(
        default=True, doc="Toggle live preview on/off", precedence=-1
    )

    color_transform_name = param.String("", doc="Name for the color transform to add")
    apply_colormap = param.Action(
        lambda self: self._action_apply_colormap(), label="Update color transform"
    )

    # show_thds = param.Boolean(True, label='Show Thresholds', doc='Toggle to show/hide threshold lines.')
    values = param.List(default=[], precedence=-1)
    colors = param.List(default=[], precedence=-1)

    def __init__(self, parent, **param):
        super(ColorTransformControl, self).__init__(
            parent, _excluded=["otsu_thd", "num_colors"], **param
        )

        # https://discourse.holoviz.org/t/based-on-a-select-widget-update-a-second-select-widget-then-how-to-link-the-latter-to-a-reactive-plot/917/8
        # update to proplot cmaps?
        self._bounds = True  # set to False to disable updating bounds on thresholds
        cc_cmaps = sorted(colorcet.cm.keys())
        mpl_cmaps = sorted(
            set(plt.colormaps()) - set("cet_" + cmap for cmap in cc_cmaps)
        )

        self._pyhdx_cmaps = {}  # Dict of pyhdx default colormaps
        f_cmap, f_norm = CMAP_NORM_DEFAULTS["foldedness"]
        self._user_cmaps = {"lily_blue": f_cmap}
        cmap_opts = [opt for opt in self.opts.values() if isinstance(opt, CmapOpts)]
        self.quantity_mapping = {}  # quantity: (cmap, norm)
        for opt in cmap_opts:
            cmap, norm = opt.cmap, opt.norm_scaled
            self._pyhdx_cmaps[cmap.name] = cmap
            field = {"dG": "dG"}.get(opt.field, opt.field)  # no longer needed
            self.quantity_mapping[field] = (cmap, norm)

        self.cmap_options = {
            "matplotlib": mpl_cmaps,  # list or dicts
            "colorcet": cc_cmaps,
            "pyhdx_default": self._pyhdx_cmaps,
            "user_defined": self._user_cmaps,
        }

        self._update_num_colors()
        self._update_num_values()
        self._update_library()

        quantity_options = [
            opt.name for opt in self.opts.values() if isinstance(opt, CmapOpts)
        ]
        self.param["quantity"].objects = quantity_options
        if self.quantity is None:
            self.quantity = quantity_options[0]

        self.update_box()

    @property
    def own_widget_names(self):
        """returns a list of names of widgets in self.widgets to be laid out in controller card"""

        initial_widgets = []
        for name in self.param:
            precedence = self.param[name].precedence
            if (precedence is None or precedence > 0) and name not in self._excluded + [
                "name"
            ]:
                initial_widgets.append(name)

        # todo control color / value fields with param.add_parameter function
        widget_names = initial_widgets + [f"value_{i}" for i in range(len(self.values))]
        if self.mode != "Colormap":
            widget_names += [f"color_{i}" for i in range(len(self.colors))]
        return widget_names

    def make_dict(self):
        return self.generate_widgets(
            num_colors=pn.widgets.IntInput,
            current_color_transform=pn.widgets.StaticText,
        )

    def get_selected_data(self):
        # todo rfu residues peptides should be expanded one more dim to add quantity column level on top
        if self.quantity is None:
            self.parent.logger.info(
                "No quantity selected to derive colormap thresholds from"
            )
            return None

        # get the table key corresponding to the selected cmap quantity
        table_key = next(
            iter((k for k in self.src.tables.keys() if k.startswith(self.quantity)))
        )

        table = self.src.get_table(table_key)
        if table is None:
            self.parent.logger.info(
                f"Table corresponding to {self.quantity!r} is empty"
            )
            return None

        field = self.opts[self.quantity].field
        if isinstance(table.columns, pd.MultiIndex):
            df = table.xs(field, level=-1, axis=1)
        else:
            df = table[field]

        opt = self.opts[self.quantity]

        return df * opt.sclf

    def get_values(self):
        """return numpy array with only the values from selected dataframe, nan omitted"""

        array = self.get_selected_data().to_numpy().flatten()
        values = array[~np.isnan(array)]

        return values

    def _action_otsu(self):
        if self.num_colors <= 1:
            return
        values = self.get_values()  # todo check for no values
        if not values.size:
            return

        # func = np.log if self.log_space else lambda x: x  # this can have NaN when in log space
        func = lambda x: x
        thds = threshold_multiotsu(func(values), classes=self.num_colors)
        widgets = [
            widget for name, widget in self.widgets.items() if name.startswith("value")
        ]
        for thd, widget in zip(thds[::-1], widgets):  # Values from high to low
            widget.start = None
            widget.end = None
            widget.value = thd  # np.exp(thd) if self.log_space else thd
        self._update_bounds()

    def _action_linear(self):
        i = 1 if self.mode == "Discrete" else 0
        values = self.get_values()
        if not values.size:
            return

        # if self.log_space:
        #     thds = np.logspace(np.log(np.min(values)), np.log(np.max(values)),
        #                        num=self.num_colors + i, endpoint=True, base=np.e)
        # else:
        thds = np.linspace(
            np.min(values), np.max(values), num=self.num_colors + i, endpoint=True
        )

        widgets = [
            widget for name, widget in self.widgets.items() if name.startswith("value")
        ]
        for thd, widget in zip(thds[i : self.num_colors][::-1], widgets):
            # Remove bounds, set values, update bounds
            widget.start = None
            widget.end = None
            widget.value = thd
        self._update_bounds()

    def _action_apply_colormap(self):
        if self.quantity is None:
            return

        cmap, norm = self.get_cmap_and_norm()
        if cmap and norm:
            cmap.name = self.color_transform_name
            # with # aggregate and execute at once: #            with param.parameterized.discard_events(opt): ??
            opt = self.opts[self.quantity]
            opt.cmap = cmap
            opt.norm_scaled = norm  # perhaps setter for norm such that externally it behaves as a rescaled thingy?

            self.quantity_mapping[self.quantity] = (cmap, norm)
            self._user_cmaps[cmap.name] = cmap

    @param.depends("colormap", "values", "colors", watch=True)
    def _preview_updated(self):
        if self.live_preview:
            pass
            # self._action_apply_colormap()

    @param.depends("quantity", watch=True)
    def _quantity_updated(self):
        cmap, norm = self.quantity_mapping[self.quantity]

        # with .. # todo accumulate events?

        preview = self.live_preview

        self.live_preview = False
        self.mode = "Colormap"

        lib = (
            "pyhdx_default" if cmap.name in self._pyhdx_cmaps.keys() else "user_defined"
        )
        self.library = lib
        self.no_coverage = to_hex(cmap.get_bad(), keep_alpha=False)
        self.colormap = cmap.name
        self.color_transform_name = cmap.name
        self.current_color_transform = cmap.name

        thds = [norm.vmax, norm.vmin]
        widgets = [
            widget for name, widget in self.widgets.items() if name.startswith("value")
        ]
        self._bounds = False  # todo decorator
        for i, (thd, widget) in enumerate(zip(thds, widgets)):
            # Remove bounds, set values, update bounds
            widget.start = None
            widget.end = None
            widget.value = thd
            self.values[i] = thd

        self.param.trigger("values")
        self._bounds = True
        self._update_bounds()
        self.live_preview = preview

    def get_cmap_and_norm(self):
        norm_klass = mpl.colors.Normalize

        # if not self.log_space else mpl.colors.LogNorm
        # if self.colormap_name in self.cmaps['pyhdx_default']: # todo change
        #     self.parent.logger.info(f"Colormap name {self.colormap_name} already exists")
        #     return None, None

        if len(self.values) < 1:
            return None, None

        if self.mode == "Discrete":
            if len(self.values) != len(self.colors) - 1:
                return None, None
            cmap = mpl.colors.ListedColormap(self.colors[::-1])
            values = self.get_values()
            thds = sorted([values.min()] + self.values + [values.max()])

            norm = mpl.colors.BoundaryNorm(
                thds, self.num_colors, extend="neither"
            )  # todo refactor values to thd_values

        elif self.mode == "Continuous":
            norm = norm_klass(
                vmin=np.min(self.values), vmax=np.max(self.values), clip=True
            )
            positions = norm(self.values[::-1])
            cmap = mpl.colors.LinearSegmentedColormap.from_list(
                "custom_cmap", list(zip(positions, self.colors))
            )

        elif self.mode == "Colormap":
            norm = norm_klass(
                vmin=np.min(self.values), vmax=np.max(self.values), clip=True
            )
            if self.library == "matplotlib":
                cmap = mpl.cm.get_cmap(self.colormap)
            elif self.library == "colorcet":
                cmap = getattr(colorcet, "m_" + self.colormap)
            elif self.library == "pyhdx_default":
                cmap = self._pyhdx_cmaps[self.colormap]
            elif self.library == "user_defined":
                try:
                    cmap = self._user_cmaps[self.colormap]
                except KeyError:
                    return None, None

        cmap.name = self.color_transform_name
        self.current_color_transform = self.color_transform_name
        cmap.set_bad(self.no_coverage)

        return cmap, norm

    @param.depends("library", watch=True)
    def _update_library(self):
        collection = self.cmap_options[self.library]
        options = (
            collection if isinstance(collection, list) else list(collection.keys())
        )
        self.param["colormap"].objects = options
        if (
            self.colormap is None or self.colormap not in options
        ) and options:  # todo how can it not be in options?
            self.colormap = options[0]

    @param.depends("mode", watch=True)
    def _mode_updated(self):
        if self.mode == "Discrete":
            self._excluded = ["library", "colormap"]
        #        self.num_colors = max(3, self.num_colors)
        #        self.param['num_colors'].bounds = (3, None)
        elif self.mode == "Continuous":
            self._excluded = ["library", "colormap", "otsu_thd"]
        #      self.param['num_colors'].bounds = (2, None)
        elif self.mode == "Colormap":
            self._excluded = ["otsu_thd", "num_colors"]
            self.num_colors = 2

        # todo adjust add/ remove color widgets methods
        self.param.trigger("num_colors")
        self.update_box()

    @param.depends("num_colors", watch=True)
    def _update_num_colors(self):
        while len(self.colors) != self.num_colors:
            if len(self.colors) > self.num_colors:
                self._remove_color()
            elif len(self.colors) < self.num_colors:
                self._add_color()
        self.param.trigger("colors")

    @param.depends("num_colors", watch=True)
    def _update_num_values(self):
        diff = 1 if self.mode == "Discrete" else 0
        while len(self.values) != self.num_colors - diff:
            if len(self.values) > self.num_colors - diff:
                self._remove_value()
            elif len(self.values) < self.num_colors - diff:
                self._add_value()

        if self.num_colors >= 5:
            self.widgets["otsu_thd"].disabled = True
        else:
            self.widgets["otsu_thd"].disabled = False

        self._update_bounds()
        self.param.trigger("values")
        self.update_box()

    def _add_value(self):
        # value widgets are ordered in decreasing order, ergo next value widget
        # starts with default value of previous value -1
        try:
            first_value = self.values[-1]
        except IndexError:
            first_value = 0

        default = float(first_value - 1)
        self.values.append(default)

        name = f"Threshold {len(self.values)}"
        key = f"value_{len(self.values) - 1}"  # values already populated, first name starts at 1
        widget = pn.widgets.FloatInput(name=name, value=default)
        self.widgets[key] = widget
        widget.param.watch(self._value_event, ["value"])

    def _remove_value(self):
        key = f"value_{len(self.values) - 1}"
        widget = self.widgets.pop(key)
        self.values.pop()

        [widget.param.unwatch(watcher) for watcher in widget.param._watchers]
        del widget

    def _add_color(self):
        try:
            default = DEFAULT_CLASS_COLORS[len(self.colors)]
        except IndexError:
            default = "#" + "".join(np.random.choice(list("0123456789abcdef"), 6))

        self.colors.append(default)

        key = f"color_{len(self.colors) - 1}"
        widget = pn.widgets.ColorPicker(value=default)

        self.widgets[key] = widget

        widget.param.watch(self._color_event, ["value"])

    def _remove_color(self):
        key = f"color_{len(self.colors) - 1}"
        widget = self.widgets.pop(key)
        self.colors.pop()
        [widget.param.unwatch(watcher) for watcher in widget.param._watchers]
        del widget

    def _color_event(self, *events):
        for event in events:
            idx = list(self.widgets.values()).index(event.obj)
            key = list(self.widgets.keys())[idx]
            widget_index = int(key.split("_")[1])
            # idx = list(self.colors_widgets).index(event.obj)
            self.colors[widget_index] = event.new

        self.param.trigger("colors")

        # todo param trigger colors????

    def _value_event(self, *events):
        """triggers when a single value gets changed"""
        for event in events:
            idx = list(self.widgets.values()).index(event.obj)
            key = list(self.widgets.keys())[idx]
            widget_index = int(key.split("_")[1])
            self.values[widget_index] = event.new

        self._update_bounds()
        self.param.trigger("values")

    def _update_bounds(self):
        # for i, widget in enumerate(self.values_widgets.values()):
        if (
            not self._bounds
        ):  # temporary fix to turn on/off bounds (perhaps should be a decorator)
            return
        for i in range(len(self.values)):
            widget = self.widgets[f"value_{i}"]
            if i > 0:
                key = f"value_{i-1}"
                prev_value = float(self.widgets[key].value)
                widget.end = np.nextafter(prev_value, prev_value - 1)
            else:
                widget.end = None

            if i < len(self.values) - 1:
                key = f"value_{i+1}"
                next_value = float(self.widgets[key].value)
                widget.start = np.nextafter(next_value, next_value + 1)
            else:
                widget.start = None


class ProteinControl(PyHDXControlPanel):

    _type = "protein"

    header = "Protein Control"

    input_mode = param.Selector(
        default="RCSB PDB Download",
        doc="Method of protein structure input",
        objects=["RCSB PDB Download", "PDB File"],
    )
    file_binary = param.Parameter(doc="Corresponds to file upload value")

    pdb_id = param.String(doc="RCSB ID of protein to download")

    load_structure = param.Action(lambda self: self._action_load_structure())

    highlight_mode = param.Selector(default="Single", objects=["Single", "Range"])

    highlight_range = param.Range(
        default=(1, 2), step=1, inclusive_bounds=[True, True],
    )
    highlight_value = param.Integer(default=1)

    highlight = param.Action(lambda self: self._action_highlight())

    clear_highlight = param.Action(lambda self: self._action_clear_highlight())

    def __init__(self, parent, **params):
        super(ProteinControl, self).__init__(
            parent, _excluded=["file_binary", "highlight_range"], **params
        )
        self.n_term, self.c_term = 1, 2
        self.src.param.watch(self._hdxm_added, "hdxm_objects")

        self.update_box()

    @property
    def _layout(self):
        return [
            ("self", self.own_widget_names),  # always use this instead of none?
            ("views.protein", "highlight_color"),
            ("transforms.protein_src", "value"),
            ("views.protein", "visual_style"),
            ("views.protein", "lighting"),
            ("views.protein", "spin"),
            ("views.protein", "reset"),
        ]

    def make_dict(self):
        return self.generate_widgets(
            file_binary=pn.widgets.FileInput(multiple=False, accept=".pdb"),
            highlight_range=pn.widgets.IntRangeSlider,
            highlight_mode=pn.widgets.RadioButtonGroup,
        )

    @param.depends("input_mode", "highlight_mode", watch=True)
    def _update_mode(self):
        _excluded = []

        if self.input_mode == "PDB File":
            _excluded.append("pdb_id")
        elif self.input_mode == "RCSB PDB Download":
            _excluded.append("file_binary")

        if self.highlight_mode == "Single":
            _excluded.append("highlight_range")
        elif self.highlight_mode == "Range":
            _excluded.append("highlight_value")

        self._excluded = _excluded
        self.update_box()

    def _action_highlight(self):
        data = {
            "color": {"r": 200, "g": 105, "b": 180},
            "focus": True,
        }

        if self.highlight_mode == "Single":
            data["residue_number"] = self.highlight_value
        elif self.highlight_mode == "Range":
            data["start_residue_number"] = self.highlight_range[0]
            data["end_residue_number"] = self.highlight_range[1]

        self.views["protein"].pdbe.highlight([data])

    def _action_clear_highlight(self):
        self.views["protein"].pdbe.clear_highlight()

    def _action_load_structure(self):
        if self.input_mode == "PDB File" and self.file_binary is not None:
            pdb_string = self.file_binary.decode()
            self.parent.sources["pdb"].add_from_string(
                pdb_string, f"local_{uuid.uuid4()}"
            )  # todo parse and extract pdb id?

        elif self.input_mode == "RCSB PDB Download":
            if len(self.pdb_id) != 4:
                self.parent.logger.info(f"Invalid RCSB pdb id: {self.pdb_id}")
                return

            self.parent.sources["pdb"].add_from_pdb(self.pdb_id)

    def _hdxm_added(self, *events):
        hdxm_object = next(reversed(self.src.hdxm_objects.values()))
        self.n_term = min(self.n_term, hdxm_object.coverage.protein.n_term)
        self.c_term = max(self.c_term, hdxm_object.coverage.protein.c_term)

        self.param["highlight_value"].bounds = (self.n_term, self.c_term)
        self.param["highlight_range"].bounds = (self.n_term, self.c_term)


class FileExportControl(PyHDXControlPanel):
    """
    <outdated docstring>
    This controller allows users to export and download datasets.

    All datasets can be exported as .txt tables.
    'Mappable' datasets (with r_number column) can be exported as .pml pymol script, which colors protein structures
    based on their 'color' column.

    """

    _type = "file_export"

    header = "File Export"

    table = param.Selector(label="Target dataset", doc="Name of the dataset to export")
    export_format = param.Selector(
        default="csv",
        objects=["csv", "pprint"],
        doc="Format of the exported tables."
        "'csv' is machine-readable, 'pprint' is human-readable format",
    )

    def __init__(self, parent, **param):
        super(FileExportControl, self).__init__(parent, **param)
        self.sources["main"].param.watch(
            self._tables_updated, ["tables", "updated"]
        )  # todo make up your mind: trigger tables or updated?

    def make_dict(self):
        widgets = self.generate_widgets()

        widgets["export_tables"] = pn.widgets.FileDownload(
            label="Download table", callback=self.table_export_callback
        )
        widgets["export_pml"] = pn.widgets.FileDownload(
            label="Download pml scripts", callback=self.pml_export_callback,
        )
        widgets["export_colors"] = pn.widgets.FileDownload(
            label="Download colors", callback=self.color_export_callback,
        )

        widget_order = [
            "table",
            "export_format",
            "export_tables",
            "export_pml",
            "export_colors",
        ]
        final_widgets = {w: widgets[w] for w in widget_order}

        return final_widgets

    def _tables_updated(self, *events):
        options = list(self.sources["main"].tables.keys())
        self.param["table"].objects = options
        if not self.table and options:
            self.table = options[0]

    @property
    def _layout(self):
        return [("self", None)]

    @param.depends("table", "export_format", watch=True)
    def _table_updated(self):

        ext = ".csv" if self.export_format == "csv" else ".txt"
        self.widgets["export_tables"].filename = self.table + ext

        qty = self.table.split("_")[0]
        cmap_opts = {
            k: opt for k, opt in self.opts.items() if isinstance(opt, CmapOpts)
        }
        if qty in cmap_opts.keys():
            self.widgets["export_pml"].disabled = False
            self.widgets["export_colors"].disabled = False
            self.widgets["export_pml"].filename = self.table + "_pml_scripts.zip"
            self.widgets["export_colors"].filename = self.table + "_colors" + ext
        else:
            self.widgets["export_pml"].disabled = True
            self.widgets["export_colors"].disabled = True

    @pn.depends("table")  # param.depends?
    def table_export_callback(self):
        if self.table:
            df = self.sources["main"].tables[self.table]
            io = dataframe_to_stringio(df, fmt=self.export_format)
            return io
        else:
            return None

    @pn.depends("table")
    def pml_export_callback(self):
        if self.table:
            # todo check if table is valid for pml conversion

            color_df = self.get_color_df()

            bio = BytesIO()
            with zipfile.ZipFile(bio, "w") as pml_zip:
                for col_name in color_df.columns:
                    name = (
                        "_".join(str(col) for col in col_name)
                        if isinstance(col_name, tuple)
                        else str(col_name)
                    )
                    colors = color_df[col_name]
                    pml_script = series_to_pymol(
                        colors
                    )  # todo refactor pd_series_to_pymol?
                    pml_zip.writestr(name + ".pml", pml_script)

            bio.seek(0)
            return bio

    def get_color_df(self):
        df = self.sources["main"].tables[self.table]
        qty = self.table.split("_")[0]
        opt = self.opts[qty]
        cmap = opt.cmap
        norm = opt.norm
        if qty == "dG":
            df = df.xs("dG", level=-1, axis=1)
        elif qty == "ddG":
            df = df.xs("ddG", level=-1, axis=1)

        color_df = apply_cmap(df, cmap, norm)

        return color_df

    @pn.depends("table")
    def color_export_callback(self):
        if self.table:
            df = self.get_color_df()
            io = dataframe_to_stringio(df, fmt=self.export_format)
            return io
        else:
            return None


class FigureExportControl(PyHDXControlPanel):

    _type = "figure_export"

    header = "Figure Export"

    figure = param.Selector(
        default="scatter", objects=["scatter", "linear_bars", "rainbowclouds"]
    )

    table = param.Selector(doc="which table data to use for figure")

    figure_selection = param.Selector(
        label="Selection", doc="for scatter / rainbowclouds, which fit to use"
    )

    groupby = param.Selector(doc="for linear bars, how to group the bars")

    reference = param.Selector(allow_None=True)

    figure_format = param.Selector(default="png", objects=["png", "pdf", "svg", "eps"])

    ncols = param.Integer(
        default=2,
        label="Number of columns",
        bounds=(1, 4),
        doc="Number of columns in subfigure",
    )

    aspect = param.Number(
        default=3.0, label="Aspect ratio", doc="Subfigure aspect ratio"
    )

    width = param.Number(
        default=cfg.getfloat("plotting", "page_width"),
        label="Figure width (mm)",
        bounds=(50, None),
        doc="""Width of the output figure""",
    )

    def __init__(self, parent, **param):
        super(FigureExportControl, self).__init__(parent, **param)
        self.sources["main"].param.watch(self._figure_updated, ["updated"])
        self._figure_updated()

    def make_dict(self):
        widgets = self.generate_widgets()

        widgets["export_figure"] = pn.widgets.FileDownload(
            label="Download figure", callback=self.figure_export_callback,
        )

        widget_order = [
            "figure",
            "table",
            "figure_selection",
            "groupby",
            "reference",
            "figure_format",
            "ncols",
            "aspect",
            "width",
            "export_figure",
        ]
        final_widgets = {w: widgets[w] for w in widget_order}

        return final_widgets

    @pn.depends("figure", watch=True)
    def _figure_updated(self, *events):
        # generalize more when other plot options are introduced
        # this needs a cross section filter probably
        if not self.figure:
            self.widgets["export_figure"].disabled = True
            return

        if self.figure == "linear_bars":
            table_options = {"dG_fits", "rfu_residues"}
        else:
            table_options = {"dG_fits"}

        options = list(table_options & self.src.tables.keys())
        self.param["table"].objects = options
        if (not self.table and options) or (options and self.table not in options):
            self.table = options[0]

        if len(options) == 0:
            self.widgets["export_figure"].disabled = True
            return

        self.widgets["export_figure"].disabled = False
        if self.figure == "rainbowclouds":
            df = self.plot_data
            options = list(df.columns.unique(level=0))
            self.param["figure_selection"].objects = options
            if (
                not self.figure_selection and options
            ):  # perhaps also check if the current selection is still in options
                self.figure_selection = options[0]

            self._update_reference()

            self.aspect = cfg.getfloat("plotting", f"{self.figure}_aspect")
            self._excluded = ["groupby", "ncols"]

        elif self.figure in ["linear_bars"]:
            # _table_updated?

            self.aspect = cfg.getfloat("plotting", f"{self.figure}_aspect")
            self._excluded = ["ncols", "figure_selection"]

        elif self.figure == "scatter":

            # move to function
            df = self.plot_data
            options = list(df.columns.unique(level=0))
            self.param["figure_selection"].objects = options
            if not self.figure_selection and options:
                self.figure_selection = options[0]

            self.aspect = cfg.getfloat("plotting", "dG_aspect")
            self._excluded = ["groupby"]
        else:
            raise ValueError("how did you manage to get here?")
            # self.aspect = cfg.getfloat('plotting', f'{self.figure}_aspect')
            # self._excluded = ['ncols']

        self.update_box()

    @property
    def plot_data(self):
        df = self.sources["main"].tables[self.table]
        return df

    @pn.depends("table", watch=True)
    def _update_groupby(self):
        # update groupby options
        options = self.plot_data.columns.names[:2]
        self.param["groupby"].objects = options
        # if self.groupby not in options:
        #     self.groupby = None
        # i guess its not set to None if its no longer in options?
        # todo investigate
        if (not self.groupby and options) or (options and self.groupby not in options):
            self.groupby = options[0]

    @pn.depends("groupby", watch=True)
    def _update_reference(self):
        # update reference options
        if self.figure == "linear_bars":
            df = self.plot_data
            groupby_index = self.plot_data.columns.names.index(self.groupby)
            barsby_index = 1 - groupby_index
            options = list(self.plot_data.columns.unique(level=barsby_index))
        else:
            options = list(self.plot_data.columns.unique(level=1))

        self.param["reference"].objects = [None] + options

    @pn.depends("figure_selection", watch=True)
    def _figure_selection_updated(self):  # selection is usually Fit ID
        df = self.sources["main"].tables["dG_fits"][self.figure_selection]
        options = list(df.columns.unique(level=0))
        self.param["reference"].objects = [None] + options
        if not self.reference and options:
            self.reference = None

    @pn.depends(
        "figure",
        "table",
        "reference",
        "groupby",
        "figure_selection",
        "figure_format",
        watch=True,
    )
    def _figure_filename_updated(self):
        if self.table == "dG_fits":
            qty = "dG" if self.reference is None else "ddG"
        elif self.table == "rfu_residues":
            qty = "rfu" if self.reference is None else "drfu"

        if self.figure == "linear_bars":
            extra = f"by_{self.groupby}"
        else:
            extra = self.figure_selection

        fname = f"{self.figure}_{qty}_{extra}.{self.figure_format}"

        self.widgets["export_figure"].filename = fname

    @pn.depends("figure")
    def figure_export_callback(self):
        self.widgets["export_figure"].loading = True

        if self.figure == "scatter":
            sub_df = self.plot_data[self.figure_selection]
            if self.reference is None:
                opts = self.opts["dG"]
                fig, axes, cbars = dG_scatter_figure(
                    sub_df, cmap=opts.cmap, norm=opts.norm, **self.figure_kwargs
                )

            else:
                opts = self.opts["ddG"]
                fig, axes, cbar = ddG_scatter_figure(
                    sub_df,
                    reference=self.reference,
                    cmap=opts.cmap,
                    norm=opts.norm,
                    **self.figure_kwargs,
                )
        elif self.figure == "linear_bars":
            if self.table == "dG_fits":
                opts = self.opts["ddG"] if self.reference else self.opts["dG"]
                field = "dG"
            elif self.table == "rfu_residues":
                opts = (
                    self.opts["drfu"] if self.reference else self.opts["rfu"]
                )  # TODO update to drfu
                field = "rfu"

            fig, axes, cbar = linear_bars_figure(
                self.plot_data,
                groupby=self.groupby,
                reference=self.reference,
                cmap=opts.cmap,
                norm=opts.norm,
                field=field,
            )
        elif self.figure == "rainbowclouds":
            sub_df = self.plot_data[self.figure_selection]
            opts = self.opts["ddG"] if self.reference else self.opts["dG"]
            fig, axes, cbar = rainbowclouds_figure(
                sub_df, reference=self.reference, cmap=opts.cmap, norm=opts.norm
            )

        bio = BytesIO()
        fig.savefig(bio, format=self.figure_format)
        bio.seek(0)

        self.widgets["export_figure"].loading = False

        return bio

    @property
    def figure_kwargs(self):
        kwargs = {
            "width": self.width,
            "aspect": self.aspect,
        }
        if self.figure == "scatter":
            kwargs["ncols"] = self.ncols
        return kwargs


class SessionManagerControl(PyHDXControlPanel):
    _type = "session_manager"

    header = "Session Manager"

    session_file = param.Parameter()

    load_session = param.Action(lambda self: self._load_session())

    reset_session = param.Action(lambda self: self._reset_session())

    def make_dict(self):
        widgets = self.generate_widgets(session_file=pn.widgets.FileInput)

        widgets["export_session"] = pn.widgets.FileDownload(
            label="Export session",
            callback=self.export_session_callback,
            filename="PyHDX_session.zip",
        )

        names = ["session_file", "load_session", "export_session", "reset_session"]
        widgets = {name: widgets[name] for name in names}

        return widgets

    def export_session_callback(self):
        dt = datetime.today().strftime("%Y%m%d_%H%M")
        self.widgets["export_session"].filename = f"{dt}_PyHDX_session.zip"
        bio = BytesIO()
        with zipfile.ZipFile(bio, "w") as session_zip:
            for name, table in self.sources["main"].tables.items():
                sio = dataframe_to_stringio(table)
                session_zip.writestr(name + ".csv", sio.getvalue())

        bio.seek(0)
        return bio

    def _load_session(self):
        if self.session_file is None:
            self.parent.logger.info("No session file selected")
            return None

        self.widgets["load_session"].loading = True
        if sys.getsizeof(self.session_file) > 5.0e8:
            self.parent.logger.info("Uploaded file is too large, maximum is 500 MB")
            return None

        bio = BytesIO(self.session_file)

        # todo: put pdb file in a source (PDBFileSource?)
        # write / read it to/from zipfile
        # make view read it
        session_zip = zipfile.ZipFile(bio)
        session_zip.printdir()
        names = set(session_zip.namelist())
        accepted_names = {
            "rfu_residues.csv",
            "rates.csv",
            "peptides.csv",
            "dG_fits.csv",
            "ddG_comparison.csv",
            "d_calc.csv",
            "loss.csv",
            "peptide_mse.csv",
        }

        self._reset()
        src = self.sources["main"]
        tables = names & accepted_names
        for name in tables:
            bio = BytesIO(session_zip.read(name))
            df = csv_to_dataframe(bio)
            df.columns = fix_multiindex_dtypes(df.columns)

            src.tables[name.split(".")[0]] = df

        src.param.trigger("tables")  # todo do not trigger tables?
        src.updated = True

        self.parent.logger.info(
            f"Successfully loaded PyHDX session file: {self.widgets['session_file'].filename}"
        )
        self.parent.logger.info(f"Containing the following tables: {', '.join(tables)}")

        self.widgets["load_session"].loading = False

    def _reset_session(self):
        self._reset()
        self.sources["main"].updated = True

        # todo for ctrl in cotrollers ctrl.reset()?

    def _reset(self):
        src = self.sources["main"]
        with param.parameterized.discard_events(src):
            src.hdxm_objects = {}
            src.rate_results = {}
            src.dG_fits = {}

        src.tables = {}  # are there any dependies on this?


class GraphControl(PyHDXControlPanel):
    _type = "graph"

    header = "Graph Control"

    def __init__(self, parent, **params):
        super(GraphControl, self).__init__(parent, **params)

        self.widget_transforms = [
            name
            for name, trs in self.transforms.items()
            if isinstance(trs, CrossSectionTransform)
        ]
        self._watchers = {}

        transforms = [self.transforms[f] for f in self.widget_transforms]
        for t in transforms:
            t.param.watch(self._transforms_redrawn, "redrawn")

    def _transforms_redrawn(self, *events):
        # todo they may all/multiple redraw at the same time so this gets called multiple times
        widgets = self.make_dict()
        transform_widgets = self._link_widgets()

        self.widgets = {**widgets, **transform_widgets}
        self.update_box()

    def _link_widgets(self):
        trs_dict = {
            flt: self.transforms[flt].widgets.keys() for flt in self.widget_transforms
        }
        grouped = {}
        for k, v in trs_dict.items():
            for vi in v:
                if vi in grouped.keys():
                    grouped[vi].append(k)
                else:
                    grouped[vi] = [k]

        output_widgets = {}

        for widget_name, trs in grouped.items():
            if len(trs) > 1:  # If there are multiple in the same group, link them
                master_widget = self.transforms[trs[0]].widgets[widget_name]
                client_widgets = [
                    self.transforms[t].widgets[widget_name] for t in trs[1:]
                ]

                for watcher in self._watchers.get(widget_name, []):
                    pass
                    # TODO check these warnings
                    # master_widget.param.unwatch(watcher)
                self._watchers[widget_name] = []

                for i, client in enumerate(client_widgets):
                    watcher = master_widget.link(client, value="value")
                    self._watchers[widget_name].append(watcher)

                output_widgets[widget_name] = master_widget
            else:
                output_widgets[widget_name] = self.transforms[trs[0]].widgets[
                    widget_name
                ]
        return output_widgets
