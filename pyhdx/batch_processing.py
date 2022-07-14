from __future__ import annotations

import warnings
from io import StringIO
from functools import reduce
from pathlib import Path
import os
import re
from typing import Union, Any, Optional

from pyhdx import TorchFitResult
from pyhdx.models import PeptideMasterTable, HDXMeasurement, HDXMeasurementSet
from pyhdx.fileIO import read_dynamx, csv_to_dataframe, save_fitresult
from pyhdx.fitting import (
    fit_rates_half_time_interpolate,
    fit_rates_weighted_average,
    fit_gibbs_global,
    fit_gibbs_global_batch,
    RatesFitResult,
    GenericFitResult,
)
from pyhdx.support import gen_subclasses

import param
import pandas as pd
import numpy as np
import yaml

time_factors = {"s": 1, "m": 60.0, "min": 60.0, "h": 3600, "d": 86400}
temperature_offsets = {"c": 273.15, "celsius": 273.15, "k": 0, "kelvin": 0}


# todo add data filters in state spec?
# todo add proline, n_term options
class StateParser(object):
    """

    Args:
        state_spec: Dictionary with HDX-MS state specification.
        data_src: Optional data source with input data files. If not specified, current
            directory is used. Otherwise, either a data source path can be specified or
            data can be given as a dictionary, where keys are filenames and values are
            :class:`~io.StringIO` with file contents.
        data_filters: Optional list of data filters to apply to the data. The filters are
            applied after the FD controls are applied and should be implemented as functions
            which take input data as a `:class:`~pandas.DataFrame` and returns the filtered
            `:class:`~pandas.DataFrame`

    """

    def __init__(
        self,
        state_spec: dict,
        data_src: Union[os.PathLike[str], str, dict[str, StringIO], None],
        data_filters: list = None,
    ) -> None:

        self.state_spec = state_spec
        data_src = data_src or "."
        if isinstance(data_src, (os.PathLike, str)):
            self.data_src = Path(data_src)
        elif isinstance(data_src, dict):
            self.data_src = data_src
        else:
            raise TypeError(
                f"Invalid data type {type(data_src)!r}, must be path or dict"
            )

        self.data_filters = data_filters or []

    def load_data(self, *filenames: os.PathLike, reader="dynamx") -> pd.DataFrame:
        if reader == "dynamx":
            read_func = read_dynamx
        else:
            raise NotImplementedError("Only reading of dynamx files is implemented")

        if isinstance(self.data_src, Path):
            input_files = [self.data_src / filename for filename in filenames]
            df = read_func(*input_files)
        else:
            input_stringios = [
                self.data_src[Path(filename).name] for filename in filenames
            ]
            df = read_func(*input_stringios)
            for io in input_stringios:
                io.seek(0)

        return df

    @property
    def hdxm_list(self) -> list[HDXMeasurement]:
        hdxm_list = []
        for state in self.state_spec.keys():
            hdxm = self.load_hdxm(state, name=state)
            hdxm_list.append(hdxm)

        return hdxm_list

    def load_hdxmset(self) -> HDXMeasurementSet:
        return HDXMeasurementSet(self.hdxm_list)

    def load_hdxm(self, state: str, **kwargs: Any) -> HDXMeasurement:
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

        # Use a FD control for back exchange correction
        if "FD_control" in state_dict:
            control_state = state_dict["FD_control"]["state"]
            exposure_value = state_dict["FD_control"]["exposure"]["value"]
            exposure_units = state_dict["FD_control"]["exposure"]["unit"]
            control_exposure = exposure_value * time_factors[exposure_units]

            if "ND_control" in state_dict:
                # TODO need a function for reading value /  unit blocks
                nd_control_state = state_dict["ND_control"]["state"]
                nd_exposure_value = state_dict["ND_control"]["exposure"]["value"]
                nd_exposure_units = state_dict["ND_control"]["exposure"]["unit"]
                nd_control_exposure = (
                    nd_exposure_value * time_factors[nd_exposure_units]
                )
                control_0 = (nd_control_state, nd_control_exposure)
            else:
                control_0 = None

            pmt.set_control(
                control_1=(control_state, control_exposure), control_0=control_0
            )
        # Flat back exchange percentage for all peptides
        elif "be_percent" in state_dict.keys():
            pmt.set_backexchange(state_dict["be_percent"])
        else:
            raise ValueError("No valid back exchange control method specified")

        if "temperature" in state_dict:
            temperature = state_dict["temperature"]["value"]
            try:
                t_offset = temperature_offsets[state_dict["temperature"]["unit"]]
            except KeyError:
                t_offset = temperature_offsets[
                    state_dict["temperature"]["unit"].lower()
                ]

            temperature += t_offset
        else:
            temperature = None

        pH = state_dict.get("pH", None)
        sequence = state_dict.get("sequence", "")
        c_term = state_dict.get("c_term")
        n_term = state_dict.get("n_term", 1)

        if not (c_term or sequence):
            raise ValueError("Must specify either 'c_term' or 'sequence'")

        state_data = pmt.get_state(state_dict["experiment"]["state"])
        if "exposure" in state_dict["experiment"]:
            t_unit = state_dict["experiment"]["exposure"]["unit"]
            times = (
                np.array(state_dict["experiment"]["exposure"]["values"])
                * time_factors[t_unit]
            )
        else:
            all_times = state_data["exposure"].unique()
            times = all_times[np.nonzero(all_times)]

        t_set = set(times)  # set of requested exposure times
        d_set = set(
            state_data["exposure"].unique()
        )  # set of exposure times in the data

        # Check if all requested exposure times are present
        if not t_set.issubset(d_set):
            diff = t_set - d_set
            raise ValueError(
                f"The following requested exposure times were not found in the "
                f"supplied data: {', '.join([str(e) for e in diff])}"
            )

            # Select only requested exposure times
        state_data = state_data[state_data["exposure"].isin(times)]

        for flt in self.data_filters:
            state_data = flt(state_data)

        if "name" not in kwargs:
            kwargs["name"] = state

        hdxm = HDXMeasurement(
            state_data,
            temperature=temperature,
            pH=pH,
            sequence=sequence,
            n_term=n_term,
            c_term=c_term,
            **kwargs,
        )

        return hdxm


process_functions = {
    "csv_to_dataframe": csv_to_dataframe,
    "fit_rates_half_time_interpolate": fit_rates_half_time_interpolate,
    "fit_rates_weighted_average": fit_rates_weighted_average,
    "fit_gibbs_global": fit_gibbs_global,
}


# task objects should be param
class Task(param.Parameterized):
    ...

    scheduler_address = param.String(doc="Optional scheduler adress for dask task")

    cwd = param.ClassSelector(Path, doc="Path of the current working directory")


class LoadHDMeasurementSetTask(Task):
    _type = "load_hdxm_set"

    state_file = param.String()  # = string path

    out = param.ClassSelector(HDXMeasurementSet)

    def execute(self, *args, **kwargs):
        state_spec = yaml.safe_load((self.cwd / self.state_file).read_text())
        parser = StateParser(state_spec, self.cwd, default_filters)
        hdxm_set = parser.load_hdxmset()

        self.out = hdxm_set


class EstimateRates(Task):
    _type = "estimate_rates"

    hdxm_set = param.ClassSelector(HDXMeasurementSet)

    select_state = param.String(
        doc="If set, only use this state for creating initial guesses"
    )

    out = param.ClassSelector((RatesFitResult, GenericFitResult))

    def execute(self, *args, **kwargs):
        if self.select_state:  # refactor to 'state' ?
            hdxm = self.hdxm_set.get(self.select_state)
            result = fit_rates_half_time_interpolate(hdxm)
        else:
            results = []
            for hdxm in self.hdxm_set:
                r = fit_rates_half_time_interpolate(hdxm)
                results.append(r)
            result = RatesFitResult(results)

        self.out = result


# todo allow guesses from deltaG
class ProcessGuesses(Task):
    _type = "create_guess"

    hdxm_set = param.ClassSelector(HDXMeasurementSet)

    select_state = param.String(
        doc="If set, only use this state for creating initial guesses"
    )

    rates_df = param.ClassSelector(pd.DataFrame)

    out = param.ClassSelector((pd.Series, pd.DataFrame))

    def execute(self, *args, **kwargs):
        if self.select_state:
            hdxm = self.hdxm_set.get(self.select_state)
            if self.rates_df.columns.nlevels == 2:
                rates_series = self.rates_df[(self.select_state, "rate")]
            else:
                rates_series = self.rates_df["rate"]

            guess = hdxm.guess_deltaG(rates_series)

        else:
            rates = self.rates_df.xs("rate", level=-1, axis=1)
            guess = self.hdxm_set.guess_deltaG(rates)

        self.out = guess


class FitGlobalBatch(Task):
    _type = "fit_global_batch"

    hdxm_set = param.ClassSelector(HDXMeasurementSet)

    initial_guess = param.ClassSelector(
        (pd.Series, pd.DataFrame), doc="Initial guesses for fits"
    )

    out = param.ClassSelector(TorchFitResult)

    def execute(self, *args, **kwargs):
        result = fit_gibbs_global_batch(self.hdxm_set, self.initial_guess, **kwargs)

        self.out = result


class SaveFitResult(Task):
    _type = "save_fit_result"

    fit_result = param.ClassSelector(TorchFitResult)

    output_dir = param.String()

    def execute(self, *args, **kwargs):
        save_fitresult(self.cwd / self.output_dir, self.fit_result)


class JobParser(object):
    """ """

    cwd = param.ClassSelector(Path, doc="Path of the current working directory")

    def __init__(self, job_spec: dict, cwd: Optional[Path] = None):
        self.job_spec = job_spec
        self.cwd = cwd or Path().cwd()

        self.tasks = {}
        self.task_classes = {
            cls._type: cls
            for cls in gen_subclasses(Task)
            if getattr(cls, "_type", None)
        }

    def resolve_var(self, var_string: str) -> Any:
        task_name, *attrs = var_string.split(".")

        return reduce(getattr, attrs, self.tasks[task_name])

    def execute(self) -> None:

        for task_spec in self.job_spec["steps"]:
            task_klass = self.task_classes[task_spec["task"]]
            skip = {"args", "kwargs", "task"}

            resolved_params = {}
            for par_name in task_spec.keys() - skip:
                value = task_spec[par_name]
                if isinstance(value, str):
                    m = re.findall(r"\$\((.*?)\)", value)
                    if m:
                        value = self.resolve_var(m[0])
                resolved_params[par_name] = value
            task = task_klass(cwd=self.cwd, **resolved_params)
            task.execute(*task_spec.get("args", []), **task_spec.get("kwargs", {}))

            self.tasks[task.name] = task


# todo configurable
default_filters = [lambda df: df.query("exposure > 0")]
