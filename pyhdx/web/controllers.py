from __future__ import annotations

import asyncio
import sys
import uuid
import zipfile
from io import StringIO, BytesIO
from typing import Any

import colorcet
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import panel as pn
import param
import yaml
from distributed import Client
from matplotlib.colors import Normalize, Colormap
from omegaconf import OmegaConf
from panel.io.server import async_execute
from proplot import to_hex
from scipy.constants import R
from skimage.filters import threshold_multiotsu

from pyhdx.config import cfg
from pyhdx.fileIO import csv_to_dataframe, dataframe_to_stringio
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
    fit_d_uptake,
    DUptakeFitResultSet,
)
from pyhdx.datasets import HDXDataSet, DataVault, DataFile
from pyhdx.fitting_torch import TorchFitResultSet
from pyhdx.models import (
    PeptideUptakeModel,
    HDXMeasurement,
)
from pyhdx.plot import (
    dG_scatter_figure,
    ddG_scatter_figure,
    linear_bars_figure,
    rainbowclouds_figure,
    CMAP_NORM_DEFAULTS,
)
from pyhdx.process import verify_sequence, correct_d_uptake, filter_peptides
from pyhdx.support import (
    series_to_pymol,
    apply_cmap,
    multiindex_astype,
    multiindex_set_categories,
    dataframe_intersection,
)
from pyhdx.web.base import ControlPanel, DEFAULT_CLASS_COLORS
from pyhdx.web.opts import CmapOpts
from pyhdx.web.sources import TABLE_INFO
from pyhdx.web.transforms import CrossSectionTransform
from pyhdx.web.utils import fix_multiindex_dtypes
from pyhdx.web.widgets import ASyncProgressBar, CompositeFloatSliders

from pyhdx.__version__ import __version__
import matplotlib

matplotlib.use("agg")


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

    field = param.String(default="init_value", doc="what is the value of this field during async?")

    @property
    def src(self):
        return self.sources["main"]

    def make_dict(self):
        widgets = self.generate_widgets()
        widgets["pbar"] = ASyncProgressBar()

        return widgets

    async def work_func(self):
        name = self.field  # is indeed stored locally
        async with Client(cfg.cluster.scheduler_address, asynchronous=True) as client:
            futures = []
            for i in range(10):
                duration = (i + np.random.rand()) / 3.0
                future = client.submit(blocking_function, duration)
                futures.append(future)

            await self.widgets["pbar"].run(futures)

            results = await asyncio.gather(*futures)

        result = pd.concat(results)
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

        input_ctrl: PeptideFileInputControl = self.parent.control_panels["PeptideFileInputControl"]
        print(f"{input_ctrl.measurement_name=}")
        print(input_ctrl.exp_file)
        print(f"{input_ctrl.exp_state=}")

    def _action_test(self):
        src = self.sources["metadata"]
        d = src.get("user_settings")

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


class GlobalSettingsControl(ControlPanel):
    _type = "global_settings"

    header = "Settings"

    drop_first = param.Integer(
        default=cfg.analysis.drop_first,
        bounds=(0, None),
        doc="Select the number of N-terminal residues to ignore.",
    )

    weight_exponent = param.Number(
        default=cfg.analysis.weight_exponent,
        bounds=(0, None),
        doc="Value of the exponent use for weighted averaging of RFU values",
    )

    def make_dict(self):
        widgets = self.generate_widgets()
        widgets["config_download"] = pn.widgets.FileDownload(
            label="Download config file", callback=self.config_download_callback
        )

        return widgets

    def config_download_callback(self) -> StringIO:
        # Generate and set filename
        timestamp = self.parent.session_time.strftime("%Y%m%d%H%M")
        self.widgets["config_download"].filename = f"PyHDX_config_{timestamp}.yaml"

        sio = StringIO()
        version_string = "# pyhdx configuration file " + __version__ + "\n\n"
        sio.write(version_string)

        masked_conf = OmegaConf.masked_copy(cfg.conf, cfg.conf.keys() - {"server"})
        OmegaConf.save(config=masked_conf, f=sio)
        sio.seek(0)

        return sio

    @param.depends("drop_first", watch=True)
    def _update_drop_first(self):
        cfg.analysis.drop_first = self.drop_first

    @param.depends("weight_exponent", watch=True)
    def _update_weight_exponent(self):
        cfg.analysis.weight_exponent = self.weight_exponent


class PeptideFileInputControl(PyHDXControlPanel):
    """
    This controller allows users to input .csv file (Currently only DynamX format) of 'state' peptide uptake data.
    Users can then choose how to correct for back-exchange and which 'state' and exposure times should be used for
    analysis.

    """

    _type = "peptide_file_input"

    header = "Peptide Input"

    input_mode = param.Selector(default="Manual", objects=["Manual", "Batch", "Database"])

    input_files_label = param.String("Input files:")

    input_files = param.List(doc="HDX input files. Currently only supports DynamX format")

    batch_file_label = param.String("Batch file (yaml)")

    batch_file = param.Parameter(doc="Batch file input:")

    dataset_id = param.Selector(
        label="Dataset ID", doc="Dataset ID to load from hdxms-datasets database"
    )

    nd_control = param.Boolean(
        default=False, precedence=-1, doc="Whether to allow users to input a ND control"
    )

    show_pH = param.Boolean(default=True, precedence=-1)

    show_temperature = param.Boolean(default=True, precedence=-1)

    show_d_percentage = param.Boolean(default=True, precedence=-1)

    fd_file = param.Selector(doc="File used for FD control", label="FD File")

    fd_state = param.Selector(doc="State used to normalize uptake", label="FD State")

    fd_exposure = param.Selector(doc="Exposure used to normalize uptake", label="FD Exposure")

    nd_file = param.Selector(doc="File used for ND control", label="ND File")

    nd_state = param.Selector(doc="State used to normalize uptake", label="ND State")

    nd_exposure = param.Selector(doc="Exposure used to normalize uptake", label="ND Exposure")

    exp_file = param.Selector(doc="File with experiment peptides", label="Exp File")

    exp_state = param.Selector(doc="State for selected experiment", label="Experiment State")

    exp_exposures = param.ListSelector(
        default=[],
        objects=[""],
        label="Experiment Exposures",
        doc="Selected exposure time to use",
    )

    d_percentage = param.Number(
        90.0,
        bounds=(0, 100),
        doc="Percentage of deuterium in the labelling buffer",
        label="Deuterium percentage",
    )

    temperature = param.Number(
        293.15,
        bounds=(273.15, 373.15),
        doc="Temperature of the D-labelling reaction",
        label="Temperature (K)",
    )

    pH = param.Number(
        7.5,
        bounds=(2.0, 14.0),
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

    measurement_name = param.String(doc="Label for the current HDX measurement")

    add_dataset_button = param.Action(  # -> refactor measurement
        lambda self: self._add_single_dataset_spec(),
        label="Add measurement",
        doc="Add single HDX measurement specification for loading",
    )

    hdxm_list = param.ListSelector(
        label="HDX Measurements", doc="Lists added HDX-MS measurements", constant=True
    )

    load_dataset_button = param.Action(
        lambda self: self._action_load_datasets(),
        label="Load dataset",
        doc="Parse specified HDX measurements apply back-exchange correction",
    )

    def __init__(self, parent, **params):
        excluded = ["batch_file", "batch_file_label"]
        super(PeptideFileInputControl, self).__init__(parent, _excluded=excluded, **params)
        self._update_mode()
        self.update_box()

        # Dictionary with current input files
        self.data_files: dict[str, DataFile] = {}
        # Dict with all files, keeps files after clearing input
        self.data_file_history: dict[str, DataFile] = {}

        # Dictionary of accumulated HDX state specifications:
        self.state_spec = {}
        # Dictionary of accumulated HDX data file specifications
        self.data_spec = {}

        # create database dir if it does not exist
        cfg.database_dir.mkdir(parents=True, exist_ok=True)
        self.data_vault = DataVault(cache_dir=cfg.database_dir)

        self.param["dataset_id"].objects = self.data_vault.datasets
        if self.data_vault.datasets:
            self.dataset_id = self.data_vault.datasets[0]

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
        widgets = self.generate_widgets(
            input_files_label=pn.widgets.StaticText(value=self.input_files_label),
            input_files=pn.widgets.FileInput(multiple=True, name="Input files"),
            batch_file_label=pn.widgets.StaticText(value=self.batch_file_label),
            batch_file=pn.widgets.FileInput(name="Batch yaml file", accept=".yaml"),
            be_percent=pn.widgets.FloatInput,
            pH=pn.widgets.FloatInput,
            temperature=pn.widgets.FloatInput,
            d_percentage=pn.widgets.FloatInput,
            sequence=text_area,
        )

        # Add hdx spec download button
        download = pn.widgets.FileDownload(
            label="Download HDX spec", callback=self.spec_download_callback
        )
        widgets["download_spec_button"] = download

        widget_order = [
            "input_mode",
            "input_files_label",
            "input_files",
            "batch_file_label",
            "batch_file",
            "dataset_id",
            "fd_file",
            "fd_state",
            "fd_exposure",
            "nd_file",
            "nd_state",
            "nd_exposure",
            "exp_file",
            "exp_state",
            "exp_exposures",
            "d_percentage",
            "temperature",
            "pH",
            "n_term",
            "c_term",
            "sequence",
            "measurement_name",
            "add_dataset_button",
            "download_spec_button",
            "hdxm_list",
            "load_dataset_button",
        ]

        sorted_widgets = {k: widgets[k] for k in widget_order}

        return sorted_widgets

    def spec_download_callback(self) -> StringIO:
        timestamp = self.parent.session_time.strftime("%Y%m%d%H%M")
        self.widgets["download_spec_button"].filename = f"PyHDX_hdx_spec_{timestamp}.yaml"

        sio = self.parent.hdx_spec_callback()
        return sio

    @param.depends("input_mode", watch=True)
    def _update_mode(self):
        excluded = set()

        # dictionary of widgets/params needed per input mode setting
        # inverse to find which widgets to exclude

        widget_dict = {
            "Manual": {
                "input_files_label",
                "input_files",
                "fd_file",
                "fd_state",
                "fd_exposure",
                "nd_state",
                "nd_exposure",
                "exp_file",
                "exp_state",
                "exp_exposures",
                "drop_first",
                "d_percentage",
                "pH",
                "temperature",
                "n_term",
                "c_term",
                "sequence",
                "add_dataset_button",
                "measurement_name",
                "download_spec_button",
            },
            "Batch": {"input_files_label", "input_files", "batch_file", "batch_file_label"},
            "Database": {
                "dataset_id",
            },
        }

        # widget_dict.pop(self.input_mode)
        excluded = (
            set.union(*(v for k, v in widget_dict.items() if k != self.input_mode))
            - widget_dict[self.input_mode]
        )
        #
        #
        # if self.input_mode == "Manual":
        #     excluded |= {"batch_file", "batch_file_label"}
        # elif self.input_mode == "Batch":
        #     excluded |= {
        #         "fd_file",
        #         "fd_state",
        #         "fd_exposure",
        #         "nd_state",
        #         "nd_exposure",
        #         "exp_file",
        #         "exp_state",
        #         "exp_exposures",
        #         "drop_first",
        #         "d_percentage",
        #         "pH",
        #         "temperature",
        #         "n_term",
        #         "c_term",
        #         "sequence",
        #         "add_dataset_button",
        #         "measurement_name",
        #         "download_spec_button",
        #     }
        # elif self.input_mode == "Database":
        #     excluded |= {
        #         "fd_file",
        #         "fd_state",
        #         "fd_exposure",
        #         "nd_state",
        #         "nd_exposure",
        #         "exp_file",
        #         "exp_state",
        #         "exp_exposures",
        #         "drop_first",
        #         "d_percentage",
        #         "pH",
        #         "temperature",
        #         "n_term",
        #         "c_term",
        #         "sequence",
        #         "add_dataset_button",
        #         "measurement_name",
        #         "download_spec_button",
        #     }

        # RFU mode input takes additional ND control
        if not self.nd_control:
            excluded |= {"nd_file", "nd_state", "nd_exposure"}

        # 'main' input mode takes additional HD labelling experiment parameters
        if not self.show_pH:
            excluded |= {"pH"}
        if not self.show_temperature:
            excluded |= {"temperature"}
        if not self.show_d_percentage:
            excluded |= {"d_percentage"}

        self._excluded = list(excluded)
        self.update_box()

    @param.depends("input_files", watch=True)
    def _read_files(self):
        if self.input_files:
            self.data_files = {
                name: DataFile(
                    name=name,
                    filepath_or_buffer=StringIO(byte_content.decode("UTF-8")),
                    format="DynamX",
                )
                for name, byte_content in zip(
                    self.widgets["input_files"].filename, self.input_files
                )
            }

            lens = [len(data_file.data) for data_file in self.data_files.values()]

            self.parent.logger.info(
                f'Loaded {len(self.input_files)} file{"s" if len(self.input_files) > 1 else ""} with a total '
                f"of {sum(lens)} peptides"
            )
        else:
            self.data_files = {}

        self.data_file_history |= self.data_files

        self.c_term = 0

        self._update_fd_file()
        self._update_fd_state()
        self._update_fd_exposure()

        self._update_nd_file()
        self._update_nd_state()
        self._update_nd_exposure()

        self._update_exp_file()
        self._update_exp_state()
        self._update_exp_exposure()

    def _update_fd_file(self):
        objects = list(self.data_files.keys())
        self.param["fd_file"].objects = objects
        self.fd_file = objects[0]

    @param.depends("fd_file", watch=True)
    def _update_fd_state(self):
        if self.data_files:
            states = list(self.data_files[self.fd_file].data["state"].unique())
            self.param["fd_state"].objects = states
            self.fd_state = states[0]
        else:
            self.param["fd_state"].objects = []

    @param.depends("fd_state", watch=True)
    def _update_fd_exposure(self):
        if self.data_files:
            df = self.data_files[self.fd_file].data
            # Get peptides only which belong to selected state
            fd_entries = df[df["state"] == self.fd_state]
            exposures = list(np.unique(fd_entries["exposure"]))
        else:
            exposures = []
        self.param["fd_exposure"].objects = exposures
        if exposures:
            self.fd_exposure = exposures[0]

    def _update_nd_file(self):
        objects = list(self.data_files.keys())
        self.param["nd_file"].objects = objects
        self.nd_file = objects[0]

    @param.depends("nd_file", watch=True)
    def _update_nd_state(self):
        if self.data_files:
            states = list(self.data_files[self.nd_file].data["state"].unique())
            self.param["nd_state"].objects = states
            self.nd_state = states[0]
        else:
            self.param["nd_state"].objects = []

    @param.depends("nd_state", watch=True)
    def _update_nd_exposure(self):
        if self.data_files:
            df = self.data_files[self.nd_file].data
            # Get peptides only which belong to selected state
            nd_entries = df[df["state"] == self.nd_state]
            exposures = list(np.unique(nd_entries["exposure"]))
        else:
            exposures = []
        self.param["nd_exposure"].objects = exposures
        if exposures:
            self.nd_exposure = exposures[0]

    def _update_exp_file(self):
        objects = list(self.data_files.keys())
        self.param["exp_file"].objects = objects
        self.exp_file = objects[0]

    @param.depends("exp_file", "fd_state", "fd_exposure", "nd_state", "nd_exposure", watch=True)
    def _update_exp_state(self):
        """Find the peptides which are both in the FD and ND states datasets, then from there determine the states which are present in the experiment dataset"""
        if self.exp_file not in self.data_files:
            self.param["exp_state"].objects = []
            return

        # IF self.has_nd... etc
        fd_spec = {"state": self.fd_state, "exposure": {"value": self.fd_exposure, "unit": "s"}}
        nd_spec = {"state": self.nd_state, "exposure": {"value": self.nd_exposure, "unit": "s"}}

        # Get the peptides which are in both the FD and ND states
        dataframes = [self.data_files[self.exp_file].data]
        dataframes.append(filter_peptides(self.data_files[self.fd_file].data, **fd_spec))
        dataframes.append(filter_peptides(self.data_files[self.fd_file].data, **nd_spec))

        intersected = dataframe_intersection(dataframes, by=["start", "stop"])
        states = list(np.unique(intersected[0]["state"]))

        self.param["exp_state"].objects = states

        # todo probably its best to clear all child selectors and then redo everything
        if self.exp_state in states:
            self._update_exp_exposure()

        elif states:
            self.exp_state = states[0]

    @param.depends("exp_state", watch=True)
    def _update_exp_exposure(self):
        if self.exp_file in self.data_files:
            df = self.data_files[self.exp_file].data
            exp_entries = df[df["state"] == self.exp_state]
            exposures = list(np.unique(exp_entries["exposure"]))
            exposures.sort()
        else:
            exposures = []

        self.param["exp_exposures"].objects = exposures
        self.exp_exposures = [e for e in exposures if e != 0.0]

        # Set default measurmenet name to the name of the state
        self.measurement_name = self.exp_state

        if not self.c_term and exposures:
            self.c_term = int(np.max(exp_entries["end"]))

    @property
    def hdx_spec(self) -> dict[str, Any]:
        return {"data_files": self.data_spec, "states": self.state_spec}

    # triggered from 'add measurement' button
    def _add_single_dataset_spec(self):
        """Adds the specifications of a single HDX Measurement to the `state_spec` / `data_spec` dictionaries"""
        if not self.data_files:
            self.parent.logger.info("No data loaded")
            return
        elif self.measurement_name in self.src.hdxm_objects.keys():
            self.parent.logger.info(f"Dataset name {self.measurement_name} already in use")
            return

        metadata = {}
        peptide_spec = {}

        exp_spec = {
            "state": self.exp_state,
            "exposure": {"values": self.exp_exposures, "unit": "s"},
        }

        peptide_spec["experiment"] = exp_spec

        df = self.data_files[self.exp_file].data
        peptides = filter_peptides(df, **exp_spec)
        corrected = correct_d_uptake(peptides)  # remove this step when _sequence field is removed
        exp_spec["data_file"] = self.exp_file

        try:
            verify_sequence(corrected, self.sequence, self.n_term, self.c_term)
        except ValueError as e:
            self.parent.logger.info(f"Cannot add dataset: {e}")
            return

        # Add the data file to the data spec
        if self.exp_file not in self.data_spec:
            self.data_spec[self.exp_file] = {
                "filename": self.exp_file,
                "format": "DynamX",
            }

        # Add the controls
        fd_spec = {
            "data_file": self.fd_file,
            "state": self.fd_state,
            "exposure": {"value": self.fd_exposure, "unit": "s"},
        }
        peptide_spec["FD_control"] = fd_spec
        if self.nd_control:
            nd_spec = {
                "data_file": self.nd_file,
                "state": self.nd_state,
                "exposure": {"value": self.nd_exposure, "unit": "s"},
            }
            peptide_spec["ND_control"] = nd_spec

        # Optionally add ph/temperature/d_percentage if this was input by the user
        if self.show_pH:
            metadata["pH"] = self.pH
        if self.show_temperature:
            metadata["temperature"] = {"value": self.temperature, "unit": "K"}
        if self.show_d_percentage:
            metadata["d_percentage"] = self.d_percentage

        metadata["n_term"] = self.n_term
        metadata["c_term"] = self.c_term
        if self.sequence:
            metadata["sequence"] = self.sequence

        self.state_spec[self.measurement_name] = {
            "peptides": peptide_spec,
            "metadata": metadata,
        }

        obj = self.param["hdxm_list"].objects or []
        self.param["hdxm_list"].objects = obj + [self.measurement_name]

    def _action_load_datasets(self) -> None:
        """Load all specified HDX measurements"""
        if self.input_mode == "Manual":
            data_src = self.data_file_history
            dataset = HDXDataSet(
                data_id=uuid.uuid4().hex, data_files=data_src, hdx_spec=self.hdx_spec
            )
        elif self.input_mode == "Batch":
            if self.hdxm_list:
                self.parent.logger.info("Cannot add data in batch after manually inputting data")
                return

            hdx_spec = yaml.safe_load(self.batch_file.decode("UTF-8"))

            # Convert loaded data_files to data src with correct keys
            data_src = {}
            for data_file, data_file_spec in hdx_spec["data_files"].items():
                data_src[data_file] = self.data_files[data_file_spec["filename"]]

            # store state spec for export
            self.state_spec = hdx_spec["states"]
            self.data_spec = hdx_spec["data_files"]

            dataset = HDXDataSet(
                data_id=uuid.uuid4().hex, data_files=data_src, hdx_spec=self.hdx_spec
            )
            self.param["hdxm_list"].objects = dataset.states
        elif self.input_mode == "Database":
            if self.dataset_id is None:
                return
            dataset = self.data_vault.load_dataset(self.dataset_id)
            self.param["hdxm_list"].objects = dataset.states
            self.parent.logger.info(f"Loaded dataset {dataset.data_id} from hdxms database")

            try:
                authors = ", ".join([author["name"] for author in dataset.metadata["authors"]])
                self.parent.logger.info(f"Author(s): {authors}")
            except KeyError:
                pass

            publications = dataset.metadata.get("publications", [])
            if publications:
                for pub in publications:
                    try:
                        pub_str = pub["title"]
                        if "DOI" in pub:
                            pub_str += f' ([{pub["DOI"]}](https://doi.org/{pub["DOI"]}))'
                        elif "URL" in pub:
                            pub_str += f' ([URL]({pub["URL"]}))'
                        self.parent.logger.info("Publication: " + pub_str)
                    except (KeyError, TypeError):
                        pass
        else:
            raise ValueError("Invalid input mode")

        # Disable input and changing config settings after loading data
        self.widgets["load_dataset_button"].disabled = True
        try:
            config_ctrl = self.parent.control_panels["GlobalSettingsControl"]
            config_ctrl.widgets["drop_first"].disabled = True
            config_ctrl.widgets["weight_exponent"].disabled = True
        except KeyError:
            pass

        for state in dataset.states:
            hdxm = HDXMeasurement.from_dataset(dataset, state)
            self.src.add(hdxm, state)
            self.parent.logger.info(
                f"Loaded experiment peptides state {hdxm.state} "
                f"({hdxm.Nt} timepoints, {len(hdxm.coverage)} peptides each)"
            )
            self.parent.logger.info(
                f"Average coverage: {hdxm.coverage.percent_coverage:.3}%, "
                f"Redundancy: {hdxm.coverage.redundancy:.2}"
            )


class DUptakeFitControl(PyHDXControlPanel):
    _type = "d_uptake_fit"

    header = "D-Uptake fit"

    repeats = param.Integer(default=25, bounds=(1, 100), doc="Number of fitting repeats")

    bounds = param.Boolean(default=True, doc="Toggle to use bounds [0 - 1]")

    r1 = param.Number(
        default=1,
        bounds=(0, None),
        doc="Value of the regularizer along residue axis.",
    )

    fit_name = param.String("D_uptake_fit_1", doc="Name for the fit result")

    _fit_names = param.List([], doc="List of current and future guess names", precedence=-1)

    do_fit = param.Action(
        lambda self: self._action_fit(),
        label="Do Fitting",
        doc="Start D-uptake fit",
    )

    def make_dict(self):
        widgets = self.generate_widgets(
            r1=pn.widgets.FloatInput,
            repeats=pn.widgets.IntInput,
        )

        widgets["pbar"] = ASyncProgressBar()

        return widgets

    def _action_fit(self):
        if len(self.src.hdxm_objects) == 0:
            self.parent.logger.info("No datasets loaded")
            return

        if self.fit_name in self._fit_names:
            self.parent.logger.info(f"D-uptake fit with name {self._fit_names} already in use")
            return

        self._fit_names.append(self.fit_name)
        self.parent.logger.info("Started D-uptake fit")
        self.param["do_fit"].constant = True
        self.widgets["do_fit"].loading = True

        user_dict = self.sources["metadata"].get("user_settings")
        user_dict["d_uptake_fit"][self.fit_name] = self.get_user_settings()
        async_execute(self._fit_d_uptake)

    def get_user_settings(self) -> dict:
        """
        Returns a dictionary with the current user settings.
        """
        keys = ["bounds", "r1"]
        d = {k: getattr(self, k) for k in keys}

        return d

    async def _fit_d_uptake(self):
        name = self.fit_name
        num_samples = len(self.src.hdxm_objects)
        guess = None

        self.widgets["pbar"].num_tasks = num_samples
        async with Client(cfg.cluster.scheduler_address, asynchronous=True) as client:
            futures = []
            for hdxm in self.src.hdxm_objects.values():
                future = client.submit(
                    fit_d_uptake,
                    hdxm,
                    guess,
                    self.r1,
                    self.bounds,
                    self.repeats,
                    False,
                    "worker_client",
                )
                futures.append(future)

            await self.widgets["pbar"].run(futures)
            results = await asyncio.gather(*futures)

        result_obj = DUptakeFitResultSet(list(results))
        self.src.add(result_obj, name)

        self.param["do_fit"].constant = False
        self.widgets["do_fit"].loading = False

        self.parent.logger.info(f"Finished D-uptake fit {name}")


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
    global_bounds = param.Boolean(default=False, doc="Set bounds globally across all datasets")
    lower_bound = param.Number(0.0, doc="Lower bound for association model fitting")

    upper_bound = param.Number(0.0, doc="Upper bound for association model fitting")

    guess_name = param.String(default="Guess_1", doc="Name for the initial guesses")

    _guess_names = param.List([], doc="List of current and future guess names", precedence=-1)

    do_fit1 = param.Action(
        lambda self: self._action_fit(),
        label="Calculate Guesses",
        doc="Start initial guess fitting",
        constant=True,
    )

    bounds = param.Dict({}, doc="Dictionary which stores rate fitting bounds", precedence=-1)

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

        user_dict = self.sources["metadata"].get("user_settings")
        user_dict["initial_guess"][self.guess_name] = self.get_user_settings()

        if self.fitting_model.lower() in ["association", "dissociation"]:
            loop = asyncio.get_running_loop()
            loop.create_task(self._fit_rates(self.guess_name))

        # this is practically instantaneous and does not require dask
        elif self.fitting_model == "Half-life (λ)":
            results = map(fit_rates_half_time_interpolate, self.src.hdxm_objects.values())

            result_obj = RatesFitResult(list(results))
            self.src.add(result_obj, self.guess_name)

            self.param["do_fit1"].constant = False
            self.widgets["do_fit1"].loading = False
            self.parent.logger.info(f"Finished initial guess fit {self.guess_name}")

    async def _fit_rates(self, name):
        num_samples = len(self.src.hdxm_objects)

        self.widgets["pbar"].num_tasks = num_samples
        async with Client(cfg.cluster.scheduler_address, asynchronous=True) as client:
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

    def get_user_settings(self) -> dict:
        """
        Returns a dictionary with the current user settings.
        """

        d = {"fitting_model": self.fitting_model}
        if self.fitting_model in ["association", "dissociation"]:
            d["global_bounds"] = self.global_bounds
            if self.global_bounds:
                d["bounds"] = [self.lower_bound, self.upper_bound]
            else:
                d["bounds"] = self.bounds

        return d


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

    epochs = param.Integer(EPOCHS, bounds=(1, None), doc="Maximum number of epochs (iterations.")

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

    _fit_names = param.List([], precedence=-1, doc="List of names of completed and running fits")

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
            self.parent.logger.info(f"Fit result with name {self.fit_name} already in use")
            return
        self._fit_names.append(self.fit_name)
        self.parent.logger.info("Started PyTorch fit")

        user_dict = self.sources["metadata"].get("user_settings")
        user_dict["dG_fit"][self.fit_name] = self.get_user_settings()

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
            dG_df = self.src.get_table("dG")

            if self.guess_mode == "One-to-one":
                gibbs_guess = dG_df.xs((self.initial_guess, "_dG"), level=[0, 2], axis=1)
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

        async with Client(cfg.cluster.scheduler_address, asynchronous=True) as client:
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
        async with Client(cfg.cluster.scheduler_address, asynchronous=True) as client:
            future = client.submit(fit_gibbs_global_batch, hdx_set, gibbs_guess, **self.fit_kwargs)
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

    def get_user_settings(self) -> dict:
        """
        Returns a dictionary with the current user settings.
        """

        d = {"initial_guess": self.initial_guess, "guess_mode": self.guess_mode}

        if self.guess_mode == "One-to-many":
            d["guess_state"] = self.guess_state
        d["fit_mode"] = self.fit_mode

        d.update(self.fit_kwargs)

        return d


class DifferentialControl(PyHDXControlPanel):
    _type = "diff"

    header = "Differential HDX"

    reference_state = param.Selector(doc="Which of the states to use as reference")

    comparison_name = param.String(default="comparison_1", doc="Name for the comparison table")

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
        # These are 'blind' transform and serve only to provide selection in this controller;
        # they are not coupled to any 'view'
        if "ddG_fit_select" in self.transforms:
            layout.append(("transforms.ddG_fit_select", None))
        if "dduptake_fit_select" in self.transforms:
            layout.append(("transforms.dduptake_fit_select", None))
        layout.append(("self", None))

        return layout

    # def get(self):
    #     print("remove this")
    #     df = self.transforms["ddG_fit_select"].get()
    #     return df

    def _source_updated(self, *events):
        # Triggered when hdxm objects are added
        options = self.src.names
        if len(options) >= 2:
            self.widgets["add_comparison"].disabled = False

        self.param["reference_state"].objects = options
        if self.reference_state is None and options:
            self.reference_state = options[0]

    def _action_add_comparison(self):
        current_df = self.src.get_table("drfu")
        if current_df is not None and self.comparison_name in current_df.columns.get_level_values(
            level=0
        ):
            self.parent.logger.info(f"Comparison name {self.comparison_name!r} already exists")
            return

        user_dict = self.sources["metadata"].get("user_settings")
        user_dict["differential_HDX"][self.comparison_name] = self.get_user_settings()

        # RFU only app has no dGs,
        if "ddG_fit_select" in self.transforms:
            self.add_ddG_comparison()
        if "dduptake_fit_select" in self.transforms:
            self.add_dd_uptake_comparison()
        self.add_drfu_comparison()

        self.parent.logger.info(f"Successfully added comparison set {self.comparison_name!r}")
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
        cov_test = dG_df.xs("covariance", axis=1, level=1).drop(self.reference_state, axis=1) ** 2
        cov = cov_test.add(cov_ref, axis=0).pow(0.5)
        columns = pd.MultiIndex.from_product(
            [[self.comparison_name], cov.columns, ["covariance"]], names=names
        )
        cov.columns = columns

        combined = pd.concat([ddG, cov], axis=1)

        categories = list(combined.columns.unique(level=1))
        combined.columns = multiindex_astype(combined.columns, 1, "category")
        combined.columns = multiindex_set_categories(combined.columns, 1, categories, ordered=True)

        self.src._add_table(combined, "ddG")

        # self.parent.sources['main'].param.trigger('tables')  #todo check/remove tables trigger

    def add_drfu_comparison(self):
        rfu_df = self.src.get_table("rfu")
        names = ["comparison_name", "comparison_state", "exposure", "quantity"]

        # Take rfu entries from df, to calculate drfu
        reference_rfu = rfu_df.xs(key=(self.reference_state, "rfu"), level=[0, 2], axis=1)
        test_rfu = rfu_df.drop(self.reference_state, axis=1, level=0).xs("rfu", level=2, axis=1)

        drfu = test_rfu.sub(reference_rfu, level="exposure").dropna(how="all", axis=1)

        # Expand multiindex level and set 'comparison_state' level as category
        columns = pd.MultiIndex.from_tuples(
            [(self.comparison_name, *cols, "drfu") for cols in drfu.columns],
            names=names,
        )
        drfu.columns = columns
        categories = list(drfu.columns.unique(level=1))
        drfu.columns = multiindex_astype(drfu.columns, 1, "category")
        drfu.columns = multiindex_set_categories(drfu.columns, 1, categories, ordered=True)

        reference_rfu_sd = rfu_df.xs(key=(self.reference_state, "rfu_sd"), level=[0, 2], axis=1)
        test_rfu_sd = rfu_df.drop(self.reference_state, axis=1, level=0).xs(
            "rfu_sd", level=2, axis=1
        )

        drfu_sd = ((test_rfu_sd**2).add((reference_rfu_sd**2))).pow(0.5).dropna(how="all", axis=1)

        # Expand multiindex level and set 'comparison_state' level as category
        columns = pd.MultiIndex.from_tuples(
            [(self.comparison_name, *cols, "drfu_sd") for cols in drfu_sd.columns],
            names=names,
        )
        drfu_sd.columns = columns
        categories = list(drfu_sd.columns.unique(level=1))
        drfu_sd.columns = multiindex_astype(drfu_sd.columns, 1, "category")
        drfu_sd.columns = multiindex_set_categories(drfu_sd.columns, 1, categories, ordered=True)

        combined = pd.concat([drfu, drfu_sd], axis=1).sort_index(axis=1)

        # TODO should be public
        self.src._add_table(combined, "drfu")

    def add_dd_uptake_comparison(self):
        d_uptake_df = self.transforms["dduptake_fit_select"].get()
        if d_uptake_df is None:
            return

        reference_d_uptake = d_uptake_df.xs(
            key=(self.reference_state, "d_uptake"), level=[0, 2], axis=1
        )
        test_d_uptake = d_uptake_df.drop(self.reference_state, axis=1, level=0).xs(
            "d_uptake", level=2, axis=1
        )

        dd_uptake = test_d_uptake.sub(reference_d_uptake, level="exposure").dropna(
            how="all", axis=1
        )

        names = ["comparison_name", "comparison_state", "exposure", "quantity"]
        columns = pd.MultiIndex.from_tuples(
            [(self.comparison_name, *cols, "dd_uptake") for cols in dd_uptake.columns],
            names=names,
        )

        dd_uptake.columns = fix_multiindex_dtypes(columns)

        self.src._add_table(dd_uptake, "dd_uptake")

    def get_user_settings(self) -> dict:
        """
        Returns a dictionary with the current user settings.
        """

        d = {"reference_state": self.reference_state}

        return d


class ColorTransformControl(PyHDXControlPanel):
    """
    This controller allows users classify 'mapping' datasets and assign them colors.

    Coloring can be either in discrete categories or as a continuous custom color map.
    """

    _type = "color_transform"

    header = "Color Transform"

    # todo unify name for target field (target_data set)
    # When coupling param with the same name together there should be an option to exclude this behaviour
    quantity = param.Selector(label="Target Quantity")  # todo refactor cmapopt / color transform??
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
    no_coverage = param.Color(default="#8c8c8c", doc="Color to use for regions of no coverage")

    live_preview = param.Boolean(default=True, doc="Toggle live preview on/off", precedence=-1)

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
        mpl_cmaps = sorted(set(plt.colormaps()) - set("cet_" + cmap for cmap in cc_cmaps))

        self._pyhdx_cmaps = {}  # Dict of pyhdx default colormaps
        f_cmap, f_norm = CMAP_NORM_DEFAULTS["foldedness"]
        self._user_cmaps = {"lily_blue": f_cmap}
        cmap_opts = [opt for opt in self.opts.values() if isinstance(opt, CmapOpts)]
        # quantity (column name / cmap_field): (cmap, norm)
        self.quantity_mapping: dict[str, (Colormap, Normalize)] = {}
        for opt in cmap_opts:
            cmap, norm = opt.cmap, opt.norm_scaled
            self._pyhdx_cmaps[cmap.name] = cmap
            # field = {"dG": "dG"}.get(opt.field, opt.field)  # no longer needed
            self.quantity_mapping[opt.field] = (cmap, norm)

        self.cmap_options = {
            "matplotlib": mpl_cmaps,  # list or dicts
            "colorcet": cc_cmaps,
            "pyhdx_default": self._pyhdx_cmaps,
            "user_defined": self._user_cmaps,
        }

        self._update_num_colors()
        self._update_num_values()
        self._update_library()

        # these are rfu, drfu, d_uptake, dg, ddg,
        quantity_options = [opt.field for opt in self.opts.values() if isinstance(opt, CmapOpts)]
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
            if (precedence is None or precedence > 0) and name not in self._excluded + ["name"]:
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
            self.parent.logger.info("No quantity selected to derive colormap thresholds from")
            return None

        # get the table key corresponding to the selected cmap quantity
        table_key = next(iter((k for k in self.src.tables.keys() if k.startswith(self.quantity))))

        table = self.src.get_table(table_key)
        if table is None:
            self.parent.logger.info(f"Table corresponding to {self.quantity!r} is empty")
            return None

        if isinstance(table.columns, pd.MultiIndex):
            df = table.xs(self.quantity, level=-1, axis=1)
        else:
            df = table[self.quantity]

        opt = self.opts[TABLE_INFO[self.quantity]["cmap_opt"]]

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
        widgets = [widget for name, widget in self.widgets.items() if name.startswith("value")]
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
        thds = np.linspace(np.min(values), np.max(values), num=self.num_colors + i, endpoint=True)

        widgets = [widget for name, widget in self.widgets.items() if name.startswith("value")]
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

            opt = self.opts[TABLE_INFO[self.quantity]["cmap_opt"]]
            opt.cmap = cmap
            opt.norm_scaled = norm  # perhaps setter for norm such that externally it behaves as a rescaled thingy?

            self.quantity_mapping[TABLE_INFO[self.quantity]["cmap_field"]] = (
                cmap,
                norm,
            )
            self._user_cmaps[cmap.name] = cmap

    @param.depends("colormap", "values", "colors", watch=True)
    def _preview_updated(self):
        if self.live_preview:
            pass
            # self._action_apply_colormap()

    @param.depends("quantity", watch=True)
    def _quantity_updated(self):
        cmap, norm = self.quantity_mapping[TABLE_INFO[self.quantity]["cmap_field"]]

        # with .. # todo accumulate events?

        preview = self.live_preview

        self.live_preview = False
        self.mode = "Colormap"

        lib = "pyhdx_default" if cmap.name in self._pyhdx_cmaps.keys() else "user_defined"
        self.library = lib
        self.no_coverage = to_hex(cmap.get_bad(), keep_alpha=False)
        self.colormap = cmap.name
        self.color_transform_name = cmap.name
        self.current_color_transform = cmap.name

        thds = [norm.vmax, norm.vmin]
        widgets = [widget for name, widget in self.widgets.items() if name.startswith("value")]
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
            norm = norm_klass(vmin=np.min(self.values), vmax=np.max(self.values), clip=True)
            positions = norm(self.values[::-1])
            cmap = mpl.colors.LinearSegmentedColormap.from_list(
                "custom_cmap", list(zip(positions, self.colors))
            )

        elif self.mode == "Colormap":
            norm = norm_klass(vmin=np.min(self.values), vmax=np.max(self.values), clip=True)
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
        options = collection if isinstance(collection, list) else list(collection.keys())
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
        if not self._bounds:  # temporary fix to turn on/off bounds (perhaps should be a decorator)
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
        default=(1, 2),
        step=1,
        inclusive_bounds=[True, True],
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
        self.n_term = min(self.n_term, hdxm_object.coverage.protein.index.min())
        self.c_term = max(self.c_term, hdxm_object.coverage.protein.index.max())

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
            label="Download pml scripts",
            callback=self.pml_export_callback,
        )
        widgets["export_colors"] = pn.widgets.FileDownload(
            label="Download colors",
            callback=self.color_export_callback,
        )

        widgets["divider"] = pn.layout.Divider()

        widgets["download_state_spec"] = pn.widgets.FileDownload(
            label="Download HDX spec",
            callback=self.hdx_spec_callback,
        )

        widgets["download_config"] = pn.widgets.FileDownload(
            label="Download config",
            callback=self.config_callback,
        )

        widgets["download_user_settings"] = pn.widgets.FileDownload(
            label="Download user settings",
            callback=self.user_settings_callback,
        )

        widgets["download_log"] = pn.widgets.FileDownload(
            label="Download log",
            callback=self.log_callback,
        )

        widget_order = [
            "table",
            "export_format",
            "export_tables",
            "export_pml",
            "export_colors",
            "divider",
            "download_state_spec",
            "download_config",
            "download_user_settings",
            "download_log",
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

        # currently only r_number indexed tables are in TABLE_INFO
        if self.table in TABLE_INFO:
            self.widgets["export_pml"].disabled = False
            self.widgets["export_colors"].disabled = False
            self.widgets["export_pml"].filename = self.table + "_pml_scripts.zip"
            self.widgets["export_colors"].filename = self.table + "_colors" + ext
        else:
            self.widgets["export_pml"].disabled = True
            self.widgets["export_colors"].disabled = True

    @pn.depends("table")
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
                    pml_script = series_to_pymol(colors)  # todo refactor pd_series_to_pymol?
                    pml_zip.writestr(name + ".pml", pml_script)

            bio.seek(0)
            return bio

    def get_color_df(self) -> pd.Dataframe:
        df = self.sources["main"].tables[self.table]
        opt = self.opts[TABLE_INFO[self.table]["cmap_opt"]]
        field = TABLE_INFO[self.table]["cmap_field"]
        cmap = opt.cmap
        norm = opt.norm
        df = df.xs(field, level=-1, axis=1)

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

    def hdx_spec_callback(self) -> StringIO:
        timestamp = self.parent.session_time.strftime("%Y%m%d%H%M")
        self.widgets["download_state_spec"].filename = f"PyHDX_hdx_spec_{timestamp}.yaml"

        sio = self.parent.hdx_spec_callback()
        return sio

    def config_callback(self) -> StringIO:
        timestamp = self.parent.session_time.strftime("%Y%m%d%H%M")
        self.widgets["download_config"].filename = f"PyHDX_config_{timestamp}.yaml"

        sio = self.parent.config_callback()
        return sio

    def user_settings_callback(self) -> StringIO:
        timestamp = self.parent.session_time.strftime("%Y%m%d%H%M")
        self.widgets["download_user_settings"].filename = f"PyHDX_config_{timestamp}.yaml"

        sio = self.parent.user_settings_callback()
        return sio

    def log_callback(self) -> StringIO:
        timestamp = self.parent.session_time.strftime("%Y%m%d%H%M")
        self.widgets["download_log"].filename = f"PyHDX_log_{timestamp}.txt"

        sio = self.parent.log_callback()
        return sio


class FigureExportControl(PyHDXControlPanel):
    _type = "figure_export"

    header = "Figure Export"

    figure = param.Selector(default="scatter", objects=["scatter", "linear_bars", "rainbowclouds"])

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

    aspect = param.Number(default=3.0, label="Aspect ratio", doc="Subfigure aspect ratio")

    width = param.Number(
        default=cfg.plotting.page_width,
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
            label="Download figure",
            callback=self.figure_export_callback,
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
            table_options = {"dG", "rfu"}
        else:
            table_options = {"dG"}

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

            self.aspect = cfg.plotting[f"{self.figure}_aspect"]
            self._excluded = ["groupby", "ncols"]

        elif self.figure in ["linear_bars"]:
            # _table_updated?

            self.aspect = cfg.plotting[f"{self.figure}_aspect"]
            self._excluded = ["ncols", "figure_selection"]

        # scatter plot is DG only
        elif self.figure == "scatter":
            # move to function
            df = self.plot_data
            options = list(df.columns.unique(level=0))
            self.param["figure_selection"].objects = options
            if not self.figure_selection and options:
                self.figure_selection = options[0]

            self.aspect = cfg.plotting.dG_aspect
            self._excluded = ["groupby"]
        else:
            raise ValueError("how did you manage to get here?")
            # self.aspect = conf.getfloat('plotting', f'{self.figure}_aspect')
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
        df = self.sources["main"].tables["dG"][self.figure_selection]
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
        if self.table == "dG":
            qty = "dG" if self.reference is None else "ddG"
        elif self.table == "rfu":
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

        # TODO use TABLE_INFO to get the correct cmap opt
        if self.figure == "scatter":
            sub_df = self.plot_data[self.figure_selection]
            if self.reference is None:
                opts = self.opts["dG_cmap"]
                fig, axes, cbars = dG_scatter_figure(
                    sub_df, cmap=opts.cmap, norm=opts.norm, **self.figure_kwargs
                )

            else:
                opts = self.opts["ddG_cmap"]
                fig, axes, cbar = ddG_scatter_figure(
                    sub_df,
                    reference=self.reference,
                    cmap=opts.cmap,
                    norm=opts.norm,
                    **self.figure_kwargs,
                )
        elif self.figure == "linear_bars":
            if self.table == "dG":
                opts = self.opts["ddG_cmap"] if self.reference else self.opts["dG_cmap"]
                field = "dG"
            elif self.table == "rfu":
                opts = (
                    self.opts["drfu_cmap"] if self.reference else self.opts["rfu_cmap"]
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
            opts = self.opts["ddG_cmap"] if self.reference else self.opts["dG_cmap"]
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
            "refaspect": self.aspect,
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
        self.widgets[
            "export_session"
        ].filename = f"{self.parent.session_time.strftime('%Y%m%d_%H%M')}_PyHDX_session.zip"
        bio = BytesIO()
        with zipfile.ZipFile(bio, "w") as session_zip:
            # Write tables
            for name, table in self.sources["main"].tables.items():
                sio = dataframe_to_stringio(table)
                session_zip.writestr(name + ".csv", sio.getvalue())

            # Write HDX measurement state specifications
            if sio := self.parent.hdx_spec_callback():
                session_zip.writestr("PyHDX_state_spec.yaml", sio.read())

            # Write config file
            sio = self.parent.config_callback()
            session_zip.writestr("PyHDX_config.yaml", sio.read())

            # Write user settings
            sio = self.parent.user_settings_callback()
            session_zip.writestr("PyHDX_user_settings.yaml", sio.read())

            # Write log file
            sio = self.parent.log_callback()
            session_zip.writestr("PyHDX_log.txt", sio.read())

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
            "rfu.csv",
            "drfu.csv",
            "rates.csv",
            "d_uptake.csv",
            "dd_uptake.csv",
            "peptides.csv",
            "dG.csv",
            "ddG.csv",
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

            table_name = name.split(".")[0]
            src.add_table(table_name, df)

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
            name for name, trs in self.transforms.items() if isinstance(trs, CrossSectionTransform)
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
        trs_dict = {flt: self.transforms[flt].widgets.keys() for flt in self.widget_transforms}
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
                client_widgets = [self.transforms[t].widgets[widget_name] for t in trs[1:]]

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
                output_widgets[widget_name] = self.transforms[trs[0]].widgets[widget_name]
        return output_widgets


DG_DEFAULT = [25.0] * 9
K_OPEN_DEFAULT = [2] * 9


class PeptidePropertiesControl(ControlPanel):
    """
    Control panel for properties of the peptide for simulating D-uptake
    """

    header = "Peptide controls"

    _type = "peptide"

    updated = param.Event(
        doc="Trigger for layout to listen to when widgets are updated",
        precedence=-1,
    )

    fasta_sequence = param.String(
        default="KLGPLTAGHH",
        doc="FASTA input for peptide sequence. First amino acid is truncated.",
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

    reload_btn = param.Action(
        lambda self: self._action_reload(),
        label="Reload peptide",
        doc="Reload with new peptide sequence",
    )

    dG = param.Array(
        default=np.array(DG_DEFAULT),
        doc="ΔG values of exchange (kJ/mol). List should be one shorter than input sequence",
    )

    k_open = param.Array(
        default=np.array(K_OPEN_DEFAULT),
        doc="Linderstrøm-Lang opening rates. Values are Log10, units s⁻¹",
    )

    k_close = param.Array(
        default=np.array([]),
        doc="Linderstrøm-Lang closing rates. Values are Log10, units s⁻¹",
    )

    dependent_variable = param.Selector(
        default="k_close",
        objects=["k_open", "k_close", "dG"],
        doc="Select which kinetic parameter should be fixed and derived from the other two",
    )

    def __init__(self, parent, **params) -> None:
        super().__init__(parent, **params)
        self._excluded = ["dG", "k_open", "k_close"]
        with param.edit_constant(self):
            self.k_close = self._get_k_close(self.dG, self.k_open)
        self.model = PeptideUptakeModel(list(self.fasta_sequence), self.temperature, self.pH)
        self.update_k_int_data()

        self.widgets = self.make_dict()  # this is the second trigger of make_dict
        self.update_box()
        self.update_d_uptake()

    @property
    def _layout(self):
        return [("self", self.own_widget_names), ("views.aa_uptake", "y")]

    def update_k_int_data(self):
        data_dict = {"aa": list(self.model.peptide), "k_int": self.model.k_int}
        df = pd.DataFrame(data_dict)
        self.src.add_table("k_int", df)

    def _action_reload(self):
        self.model = PeptideUptakeModel(list(self.fasta_sequence), self.temperature, self.pH)
        self.update_k_int_data()

        self.dependent_variable = "k_close"
        with param.parameterized.batch_call_watchers(self):
            self.dG = np.array(
                DG_DEFAULT[: len(self.model)] + [35.0] * (len(self.model) - len(DG_DEFAULT))
            )
            self.k_open = np.array(
                K_OPEN_DEFAULT[: len(self.model)] + [2.0] * (len(self.model) - len(K_OPEN_DEFAULT))
            )

        with param.edit_constant(self):
            self.k_close = self._get_k_close(self.dG, self.k_open)
        self.widgets = self.make_dict()
        self.updated = True
        self.update_box()
        self.update_d_uptake()

    def make_dict(self):
        dG_limits = {"start": 10, "end": 50}  # kJ/mol
        k_open_limits = {"start": -2, "end": 4}  # Log10
        k_close_limits = {
            k: self._get_k_close(dG_limits[k], k_open_limits[k]) for k in dG_limits.keys()
        }

        model = getattr(self, "model", None)
        names = model.peptide if model is not None else []

        widget_spec = dict(
            widget_type=CompositeFloatSliders,
            orientation="vertical",
            names=names,
        )

        widgets = self.generate_widgets(
            temperature=pn.widgets.FloatInput,
            dG={**widget_spec, **dG_limits},
            k_open={**widget_spec, **k_open_limits},
            k_close={**widget_spec, "disabled": True, **k_close_limits},
        )

        return widgets

    @param.depends("dG", watch=True)
    def value_updated(self):
        if self.dependent_variable == "dG":
            return
        elif self.dependent_variable == "k_open":
            self.k_open = self._get_k_open(self.dG, self.k_close)
            self.update_d_uptake()
        elif self.dependent_variable == "k_close":
            self.k_close = self._get_k_close(self.dG, self.k_open)
            self.update_d_uptake()

    @param.depends("k_open", watch=True)
    def k_open_updated(self):
        if self.dependent_variable == "k_open":
            return
        elif self.dependent_variable == "dG":
            self.dG = self._get_dG(self.k_open, self.k_close)
            self.update_d_uptake()
        elif self.dependent_variable == "k_close":
            self.k_close = self._get_k_close(self.dG, self.k_open)
            self.update_d_uptake()

    @param.depends("k_close", watch=True)
    def k_close_updated(self):
        if self.dependent_variable == "k_close":
            return
        elif self.dependent_variable == "dG":
            self.dG = self._get_dG(self.k_open, self.k_close)
            self.update_d_uptake()
        elif self.dependent_variable == "k_open":
            self.k_open = self._get_k_open(self.dG, self.k_close)
            self.update_d_uptake()

    @param.depends("dependent_variable", watch=True)
    def _fixed_quantity_updated(self):
        widget_keys = ["dG", "k_open", "k_close"]

        widget_keys.remove(self.dependent_variable)
        self.widgets[self.dependent_variable].disabled = True
        for widget_key in widget_keys:
            self.widgets[widget_key].disabled = False

    def update_d_uptake(self):
        time = np.logspace(-2, 6, num=250)

        d_uptake = self.model.eval_analytical(time, 10.0**self.k_open, 10.0**self.k_close)

        cols = [f"aa_{i}" for i in range(len(self.model))]
        idx = pd.Index(time, name="time")
        df = pd.DataFrame(d_uptake, index=idx, columns=cols)
        df["sum"] = df.sum(axis=1)

        self.src.add_table("d_uptake", df)
        self.src.updated = True

    # TODO input can also be numpy arrays
    def _get_k_open(self, dG: npt.ArrayLike, k_close: npt.ArrayLike) -> npt.ArrayLike:
        return k_close - dG * 1e3 / (np.log(10) * R * self.temperature)

    def _get_k_close(self, dG: npt.ArrayLike, k_open: npt.ArrayLike) -> npt.ArrayLike:
        return k_open + dG * 1e3 / (np.log(10) * R * self.temperature)

    def _get_dG(self, k_open: npt.ArrayLike, k_close: npt.ArrayLike) -> npt.ArrayLike:
        return np.log(10) * (k_close - k_open) * 1e-3 * R * self.temperature

    # move to base class?
    @property
    def src(self):
        return self.sources["main"]

    def update_limits(self):
        ...
        # k_close = k_open * np.exp(dG / (R * temp))

        # dG = np.log(k_close / k_open) * (8.3 * temp)


# k_open. dG, log10 kclose, kclose
# 0.01 10000.0 -0.21819477668335818 0.6050694463640514
# 0.01 50000.0 6.90902611658321 8110098.269947465
# 10000.0 10000.0 5.781805223316642 605069.4463640513
# 10000.0 50000.0 12.90902611658321 8110098269947.465
