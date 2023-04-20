from __future__ import annotations

from pathlib import Path
from typing import Optional, Dict, List, TYPE_CHECKING

import pandas as pd

from pyhdx.support import multiindex_set_categories, multiindex_astype

if TYPE_CHECKING:
    from pyhdx.web.controllers import PeptideFileInputControl

# todo merge with batch_processing
time_factors = {"s": 1, "m": 60.0, "min": 60.0, "h": 3600, "d": 86400}
temperature_offsets = {"c": 273.15, "celsius": 273.15, "k": 0, "kelvin": 0}


def get_view(
    widget,
):  # widget is viewable or widget or custom somthing with view method
    """If the widget object has a `view` attribute, return this, otherwise return the widget"""
    if hasattr(widget, "view"):
        return widget.view
    else:
        return widget


def load_state(
    file_input: PeptideFileInputControl,
    hdx_spec: Dict,
    data_dir: Path,
    states: Optional[List[str]] = None,
    names: Optional[List[str]] = None,
) -> None:
    """Loader method for `PeptideFileInputControl` from a HDX-MS state specification."""

    input_files = [data_dir / f_dict["filename"] for f_dict in hdx_spec["data_files"].values()]
    file_input.widgets["input_files"].filename = [f.name for f in input_files]
    f_bytes = [f.read_bytes() for f in input_files]
    file_input.input_files = f_bytes

    # file_input._read_files() # triggered by setting input files
    state_spec = hdx_spec["states"]
    states = states or list(state_spec.keys())
    names = names or states

    for state, name in zip(states, names):
        peptide_spec = state_spec[state]["peptides"]
        data_spec = hdx_spec["data_files"]

        file_input.fd_file = data_spec[peptide_spec["FD_control"]["data_file"]]["filename"]
        file_input.fd_state = peptide_spec["FD_control"]["state"]
        file_input.fd_exposure = (
            peptide_spec["FD_control"]["exposure"]["value"]
            * time_factors[peptide_spec["FD_control"]["exposure"]["unit"]]
        )

        if "ND_control" in peptide_spec:
            file_input.nd_file = data_spec[peptide_spec["ND_control"]["data_file"]]["filename"]
            file_input.nd_state = peptide_spec["ND_control"]["state"]
            file_input.nd_exposure = (
                peptide_spec["ND_control"]["exposure"]["value"]
                * time_factors[peptide_spec["ND_control"]["exposure"]["unit"]]
            )

        file_input.exp_file = data_spec[peptide_spec["experiment"]["data_file"]]["filename"]
        file_input.exp_state = peptide_spec["experiment"]["state"]

        # Set experiment exposures if specified, otherwise leave default selection
        try:
            exp_vals = peptide_spec["experiment"]["exposure"]["values"]
            f = time_factors[peptide_spec["experiment"]["exposure"]["unit"]]
            file_input.exp_exposures = [v * f for v in exp_vals]
        except KeyError:
            pass

        file_input.measurement_name = name

        metadata = state_spec[state]["metadata"]
        file_input.pH = metadata["pH"]

        file_input.temperature = (
            metadata["temperature"]["value"]
            + temperature_offsets[metadata["temperature"]["unit"].lower()]
        )

        file_input.d_percentage = metadata["d_percentage"]

        if "n_term" in metadata:
            file_input.n_term = metadata["n_term"]
        if "c_term" in metadata:
            file_input.c_term = metadata["c_term"]
        if "sequence" in metadata:
            file_input.sequence = metadata["sequence"]

        file_input._add_single_dataset_spec()


def fix_multiindex_dtypes(index: pd.MultiIndex) -> pd.MultiIndex:
    """Assigns correct dtypes to (column) multiindex"""
    for level, name in enumerate(index.names):
        if name in [
            "state",
            "D_uptake_fit_ID",
            "guess_ID",
            "fit_ID",
            "comparison_name",
            "comparison_state",
        ]:
            index = multiindex_astype(index, level, "category")
            categories = list(index.unique(level=level))
            index = multiindex_set_categories(index, level, categories, ordered=True)
        elif name in ["exposure"]:
            index = multiindex_astype(index, level, "float")

    return index
