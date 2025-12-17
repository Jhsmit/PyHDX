from __future__ import annotations
from typing import Any

from pyhdx import HDXMeasurement
from pyhdx.batch_processing import StateParser, time_factors, temperature_offsets


from typing import Literal, Optional, Union

import pandas as pd


def legacy_parsers(version):
    loader = LOADER_VERSIONS[version]
    parser = type(f"StateParser_{version}", (StateParser,), {"load_hdxm": loader})

    return parser


def load_hdxm_v041(self, state: str, **kwargs: Any) -> HDXMeasurement:
    """Read a single protein state to :class:`~pyhdx.models.HDXMeasurement`.

    Args:
        state: Name of the protein state to read.
        **kwargs: Additional keyword arguments passed to :class:`~pyhdx.models.HDXMeasurement`.

    Returns:
        The requested :class:`~pyhdx.models.HDXMeasurement`.

    """

    state_dict = self.state_spec[state]

    filenames = state_dict["filenames"]
    df = self.load_data(*filenames)

    pmt = PeptideMasterTable(
        df,
        drop_first=state_dict.get("drop_first", 1),
        d_percentage=state_dict["d_percentage"],
    )

    if "control" in state_dict.keys():  # Use a FD control for back exchange correction
        # todo control should be set from an external file
        control_state = state_dict["control"]["state"]
        exposure_value = state_dict["control"]["exposure"]["value"]
        exposure_units = state_dict["control"]["exposure"]["unit"]
        control_exposure = exposure_value * time_factors[exposure_units]

        pmt.set_control((control_state, control_exposure))
    elif "be_percent" in state_dict.keys():  # Flat back exchange percentage for all peptides\
        pmt.set_backexchange(state_dict["be_percent"])
    else:
        raise ValueError("No valid back exchange control method specified")

    temperature = state_dict["temperature"]["value"]
    try:
        t_offset = temperature_offsets[state_dict["temperature"]["unit"]]
    except KeyError:
        t_offset = temperature_offsets[state_dict["temperature"]["unit"].lower()]

    temperature += t_offset

    sequence = state_dict.get("sequence", "")
    c_term = state_dict.get("c_term")
    n_term = state_dict.get("n_term", 1)

    if not (c_term or sequence):
        raise ValueError("Must specify either 'c_term' or 'sequence'")

    state_data = pmt.get_state(state_dict["state"])
    for flt in self.data_filters:
        state_data = flt(state_data)

    if "name" not in kwargs:
        kwargs["name"] = state

    hdxm = HDXMeasurement(
        state_data,
        temperature=temperature,
        pH=state_dict["pH"],
        sequence=sequence,
        n_term=n_term,
        c_term=c_term,
        **kwargs,
    )

    return hdxm


TIME_FACTORS = {"s": 1, "m": 60.0, "min": 60.0, "h": 3600, "d": 86400}
TEMPERATURE_OFFSETS = {"c": 273.15, "celsius": 273.15, "k": 0.0, "kelvin": 0.0}


# overload typing to get correct return type
def convert_temperature(
    temperature_dict: dict, target_unit: str = "c"
) -> Union[float, list[float]]:
    """
    Convenience function to convert temperature values.

    Args:
        temperature_dict: Dictionary with temperature value(s) and unit.
        target_unit: Target unit for temperature. Must be "c", "k", "celsius", or "kelvin" and is
            case-insensitive.

    Returns:
        Converted temperature value(s).
    """

    src_unit = temperature_dict["unit"].lower()
    temp_offset = TEMPERATURE_OFFSETS[src_unit] - TEMPERATURE_OFFSETS[target_unit.lower()]
    if values := temperature_dict.get("values"):
        return [v + temp_offset for v in values]
    elif value := temperature_dict.get("value"):
        return value + temp_offset
    else:
        raise ValueError("Invalid temperature dictionary")


def convert_time(
    time_dict: dict, target_unit: Literal["s", "min", "h"] = "s"
) -> Union[float, list[float]]:
    """
    Convenience function to convert time values.

    Args:
        time_dict: Dictionary with time value(s) and unit.
        target_unit: Target unit for time.

    Returns:
        Converted time value(s).
    """

    src_unit = time_dict["unit"]

    time_factor = TIME_FACTORS[src_unit] / TIME_FACTORS[target_unit]
    if values := time_dict.get("values"):
        return [v * time_factor for v in values]
    elif value := time_dict.get("value"):
        return value * time_factor
    else:
        raise ValueError("Invalid time dictionary")


def filter_peptides(
    df: pd.DataFrame,
    state: Optional[str] = None,
    exposure: Optional[dict] = None,
    query: Optional[list[str]] = None,
    dropna: bool = True,
    time_unit: str = "s",
) -> pd.DataFrame:
    """
    Convenience function to filter a peptides DataFrame. .

    Args:
        df: Input dataframe.
        state: Name of protein state to select.
        exposure: Exposure value(s) to select. Exposure is given as a :obj:`dict`, with keys "value" or "values" for
            exposure value, and "unit" for the time unit.
        query: Additional queries to pass to [pandas.DataFrame.query][].
        dropna: Drop rows with `NaN` uptake entries.
        time_unit: Time unit for exposure column of supplied dataframe.

    Examples:
        Filter peptides for a specific protein state and exposure time:

        >>> d = {"state", "SecB WT apo", "exposure": {"value": 0.167, "unit": "min"}
        >>> filtered_df = filter_peptides(df, **d)

    Returns:
        Filtered dataframe.
    """

    if state is not None:
        df = df[df["state"] == state]

    if exposure is not None:
        t_val = convert_time(exposure, time_unit)  # type: ignore
        if isinstance(t_val, list):
            df = df[df["exposure"].isin(t_val)]
        else:
            df = df[df["exposure"] == t_val]

    if query:
        for q in query:
            df = df.query(q)

    if dropna:
        df = df.dropna(subset=["uptake"])

    return df.reset_index(drop=True)


LOADER_VERSIONS = {"041": load_hdxm_v041}
