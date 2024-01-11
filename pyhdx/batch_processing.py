from __future__ import annotations

import os
from dataclasses import dataclass
from functools import cached_property
from io import StringIO
from pathlib import Path
from typing import Union, Literal, Optional
import warnings

import pandas as pd

from pyhdx.config import cfg
from pyhdx.fileIO import read_dynamx
from pyhdx.models import HDXMeasurement, HDXMeasurementSet
from pyhdx.process import correct_d_uptake, apply_control


time_factors = {"s": 1, "m": 60.0, "min": 60.0, "h": 3600, "d": 86400}
temperature_offsets = {"c": 273.15, "celsius": 273.15, "k": 0, "kelvin": 0}


@dataclass(frozen=True)
class DataFile(object):
    name: str

    format: Literal["DynamX"]

    filepath_or_buffer: Union[Path, StringIO]

    def __post_init__(self):
        warnings.warn(
            "Will be removed in favour of the `hdxms-datasets` package ", DeprecationWarning
        )

    @cached_property
    def data(self) -> pd.DataFrame:
        if self.format == "DynamX":
            data = read_dynamx(self.filepath_or_buffer)
        else:
            raise ValueError(f"Invalid format {self.format!r}")

        if isinstance(self.filepath_or_buffer, StringIO):
            self.filepath_or_buffer.seek(0)

        return data


class StateParser(object):
    """

    Args:
        hdx_spec: Dictionary with HDX-MS state specification.
        data_src: Optional data source with input data files. If not specified, current
            directory is used. Otherwise, either a data source path can be specified or
            data can be given as a dictionary, where keys are filenames and values are
            :class:`~io.StringIO` with file contents.
    """

    def __init__(
        self,
        hdx_spec: dict,
        data_src: Union[os.PathLike[str], str, dict[str, DataFile], None],
        # filter_kwargs: Optional[dict[str, Any]] = None,
        # correction_kwargs: Optional[dict[str, Any]] = None,
    ) -> None:
        warnings.warn(
            "Will be removed in favour of the `hdxms-datasets` package ", DeprecationWarning
        )
        self.hdx_spec = hdx_spec
        self.data_files: dict[str, DataFile] = {}

        if isinstance(data_src, (os.PathLike, str)):
            data_src = Path(data_src) or Path(".")
            for name, spec in self.hdx_spec["data_files"].items():
                datafile = DataFile(
                    name=name,
                    filepath_or_buffer=data_src / spec["filename"],
                    **{k: v for k, v in spec.items() if k != "filename"},
                )
                self.data_files[name] = datafile

        elif isinstance(data_src, dict):
            self.data_files = data_src
        else:
            raise TypeError(f"Invalid data type {type(data_src)!r}, must be path or dict")

    def load_hdxmset(self) -> HDXMeasurementSet:
        hdxm_list = [self.load_hdxm(state) for state in self.hdx_spec["states"].keys()]
        return HDXMeasurementSet(hdxm_list)

    def load_peptides(self, state: Union[str, int], peptides: str) -> pd.DataFrame:
        state = self.states[state] if isinstance(state, int) else state
        peptide_spec = self.hdx_spec["states"][state]["peptides"][peptides]

        df = self.data_files[peptide_spec["data_file"]].data

        # filter_fields = {"state", "exposure", "query", "dropna"}
        # peptide_df = filter_peptides(
        #     df, **{k: v for k, v in peptide_spec.items() if k in filter_fields}
        # )

        filter_fields = {"state", "exposure", "query", "dropna"}
        peptide_df = batch_filter_peptides(
            df, **{k: v for k, v in peptide_spec.items() if k in filter_fields}
        )

        return peptide_df

    # -> function as monkey patch dataset parser; OR perhaps add them to internal dict of loaders ?
    def load_hdxm(self, state: Union[str, int]) -> HDXMeasurement:
        state = self.states[state] if isinstance(state, int) else state
        peptide_spec = self.hdx_spec["states"][state]["peptides"]
        metadata = self.hdx_spec["states"][state]["metadata"]

        peptides = self.load_peptides(state, "experiment")
        fd_peptides = (
            self.load_peptides(state, "FD_control") if "FD_control" in peptide_spec else None
        )
        nd_peptides = (
            self.load_peptides(state, "ND_control") if "ND_control" in peptide_spec else None
        )

        if fd_peptides is None and "be_percent" in metadata:
            peptides = correct_d_uptake(peptides, d_percentage=metadata.get("d_percentage", 100.0))
            back_exchange = metadata["be_percent"] / 100.0
            peptides["rfu"] = peptides["uptake"] / ((1 - back_exchange) * peptides["ex_residues"])
            peptides["uptake_corrected"] = peptides["uptake"] / (1 - back_exchange)
        elif isinstance(fd_peptides, pd.DataFrame):
            peptides = apply_control(peptides, fd_peptides, nd_peptides)
            peptides = correct_d_uptake(
                peptides,
                drop_first=cfg.analysis.drop_first,
                d_percentage=metadata.get("d_percentage", 100.0),
            )

        global_metadata = self.hdx_spec.get("metadata", {})
        global_metadata.update(metadata)
        hdxm = HDXMeasurement(peptides, name=state, **global_metadata)

        return hdxm

    @property
    def correction_kwargs(self):
        kwargs = {
            "drop_first": cfg.analysis.drop_first,
            "d_percentage": self.hdx_spec["metadata"].get("d_percentage", 100.0),
        }

        # todo:
        # if 'corrections' in self.hdx_spec:
        # ...

        return kwargs

    @property
    def states(self) -> list[str]:
        return list(self.hdx_spec["states"].keys())


# borrowed from hdxms-datasets
def batch_filter_peptides(
    df: pd.DataFrame,
    state: Optional[str] = None,
    exposure: Optional[dict] = None,
    query: Optional[list[str]] = None,
    dropna: bool = True,
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

    Examples:
        Filter peptides for a specific protein state and exposure time:

        >>> d = {"state", "SecB WT apo", "exposure": {"value": 0.167, "unit": "min"}
        >>> filtered_df = filter_peptides(df, **d)

    Returns:
        Filtered dataframe.
    """
    warnings.warn("Will be removed in favour of the `hdxms-datasets` package ", DeprecationWarning)

    if state is not None:
        df = df[df["state"] == state]

    if exposure is not None:
        t_val = batch_convert_time(exposure, target_unit="s")
        if isinstance(t_val, list):
            df = df[df["exposure"].isin(t_val)]
        else:
            df = df[df["exposure"] == t_val]

    if query:
        for q in query:
            df = df.query(q)

    if dropna:
        df = df.dropna(subset=["uptake"])

    return df


def batch_convert_time(
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

    warnings.warn("Will be removed in favour of the `hdxms-datasets` package ", DeprecationWarning)
    src_unit = time_dict["unit"]

    time_factor = time_factors[src_unit] / time_factors[target_unit]
    if values := time_dict.get("values"):
        return [v * time_factor for v in values]
    elif value := time_dict.get("value"):
        return value * time_factor
    else:
        raise ValueError("Invalid time dictionary")
