from pathlib import Path

import pandas as pd

from pyhdx.support import multiindex_set_categories, multiindex_astype

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


def load_state(ctrl, state_spec, data_dir, name=None):
    """
    Load a HDXMeasurement into the web interface
    Experimental use only

    Parameters
    ----------
    ctrl: :class:`~pyhdx.web.main_controllers.PyHDXController`
        Main controller to load the data into
    state_spec: :obj:`dict`
        Dictionary with measurement info according to batch processing format
    data_dir: path_like
    name: :obj:`str`, optional
        Optional name for the HDXMeasurement

    Returns
    -------
        None
    """
    # raise DeprecationWarning("This should just call the load_from_spec function on controller")

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in state_spec["filenames"]]
    else:
        input_files = [Path(p) for p in state_spec["filenames"]]

    files = [f.read_bytes() for f in input_files]

    file_input = ctrl.control_panels["PeptideFileInputControl"]
    file_input.input_files = files
    file_input.widgets["input_files"].filename = state_spec["filenames"]

    control_state = state_spec["FD_control"]["state"]
    control_exp = (
        state_spec["FD_control"]["exposure"]["value"]
        * time_factors[state_spec["FD_control"]["exposure"]["unit"]]
    )

    file_input.fd_state = control_state
    file_input.fd_exposure = control_exp
    file_input.pH = state_spec["pH"]
    file_input.temperature = (
        state_spec["temperature"]["value"]
        + temperature_offsets[state_spec["temperature"]["unit"].lower()]
    )
    file_input.d_percentage = state_spec["d_percentage"]

    file_input.exp_state = state_spec["experiment"]["state"]

    file_input.measurement_name = name or state_spec["state"]

    # file_input._action_load_datasets()


def load_state_rfu(ctrl, state_spec, data_dir, name=None):
    """Loader counterpart for RFU app. Even more experimental than `load_state`"""

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in state_spec["filenames"]]
    else:
        input_files = [Path(p) for p in state_spec["filenames"]]

    files = [f.read_bytes() for f in input_files]

    file_input = ctrl.control_panels["PeptideRFUFileInputControl"]
    file_input.input_files = files
    file_input.widgets["input_files"].filename = state_spec["filenames"]

    file_input.fd_state = state_spec["FD_control"]["state"]
    file_input.fd_exposure = (
        state_spec["FD_control"]["exposure"]["value"]
        * time_factors[state_spec["FD_control"]["exposure"]["unit"]]
    )

    file_input.nd_state = state_spec["ND_control"]["state"]
    file_input.nd_exposure = (
        state_spec["ND_control"]["exposure"]["value"]
        * time_factors[state_spec["ND_control"]["exposure"]["unit"]]
    )

    file_input.d_percentage = state_spec["d_percentage"]
    file_input.exp_state = state_spec["experiment"]["state"]
    file_input.measurement_name = name or state_spec["experiment"]["state"]


def fix_multiindex_dtypes(index: pd.MultiIndex) -> pd.MultiIndex:
    """Assigns correct dtypes to (column) multiindex"""
    if index.names[0] in ["state", "guess_ID", "fit_ID"]:
        index = multiindex_astype(index, 0, "category")
        categories = list(index.unique(level=0))
        index = multiindex_set_categories(index, 0, categories, ordered=True)

    if "exposure" in index.names:
        level = index.names.index("exposure")
        index = multiindex_astype(index, level, "float")

    return index
