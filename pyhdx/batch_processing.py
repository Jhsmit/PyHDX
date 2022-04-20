from __future__ import annotations

import warnings
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
import param
import pandas as pd
from pyhdx.support import gen_subclasses
import yaml

time_factors = {"s": 1, "m": 60.0, "min": 60.0, "h": 3600, "d": 86400}
temperature_offsets = {"c": 273.15, "celsius": 273.15, "k": 0, "kelvin": 0}


# todo add data filters in state spec?
# todo add proline, n_term options
class StateParser(object):
    "" "object used to parse yaml state input files into PyHDX HDX Measurement object"

    def __init__(
        self,
        state_spec: dict,
        data_src: Union[os.PathLike, dict, None],
        data_filters: list = None,
    ):

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

    def load_hdxmset(self) -> HDXMeasurementSet:
        """batch read the full yaml spec into a hdxmeasurementset"""
        hdxm_list = []
        for state in self.state_spec.keys():
            hdxm = self.load_hdxm(state, name=state)
            hdxm_list.append(hdxm)

        return HDXMeasurementSet(hdxm_list)

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

        if (
            "control" in state_dict.keys()
        ):  # Use a FD control for back exchange correction
            # todo control should be set from an external file
            control_state = state_dict["control"]["state"]
            exposure_value = state_dict["control"]["exposure"]["value"]
            exposure_units = state_dict["control"]["exposure"]["unit"]
            control_exposure = exposure_value * time_factors[exposure_units]

            pmt.set_control((control_state, control_exposure))
        elif (
            "be_percent" in state_dict.keys()
        ):  # Flat back exchange percentage for all peptides\
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
        n_term = state_dict.get("n_term") or 1

        if not (c_term or sequence):
            raise ValueError("Must specify either 'c_term' or 'sequence'")

        state_data = pmt.get_state(state_dict["state"])
        for flt in self.data_filters:
            state_data = flt(state_data)

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

    cwd = param.ClassSelector(Path, doc="Path of the current working directory")

    def __init__(self, job_spec: dict, cwd: Optional[os.PathLike] = None):
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


def yaml_to_hdxmset(yaml_dict, data_dir=None, **kwargs):
    """reads files according to `yaml_dict` spec from `data_dir into HDXMEasurementSet"""

    warnings.warn("yaml_to_hdxmset is deprecated, use 'StateParser'")
    hdxm_list = []
    for k, v in yaml_dict.items():
        hdxm = yaml_to_hdxm(v, data_dir=data_dir, name=k)
        hdxm_list.append(hdxm)

    return HDXMeasurementSet(hdxm_list)


# todo configurable
default_filters = [lambda df: df.query("exposure > 0")]


def yaml_to_hdxm(yaml_dict, data_dir=None, data_filters=None, **kwargs):
    # todo perhas classmethod on HDXMeasurement object?
    # merge with method in
    """
    Creates a :class:`~pyhdx.models.HDXMeasurement` object from dictionary input.

    Dictionary can be generated from .yaml format. See templates/yaml_files/SecB.yaml for format specification.

    Parameters
    ----------
    yaml_dict : :obj:`dict`
        Input dictionary specifying experimental metadata and file location to load
    data_dir : :obj:`str` or pathlib.Path object

    Returns
    -------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
        Output data object as specified by `yaml_dict`.
    """

    warnings.warn(
        "This method is deprecated in favor of StateParser", DeprecationWarning
    )

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in yaml_dict["filenames"]]
    else:
        input_files = yaml_dict["filenames"]

    data = read_dynamx(*input_files)

    pmt = PeptideMasterTable(
        data,
        drop_first=yaml_dict.get("drop_first", 1),
        d_percentage=yaml_dict["d_percentage"],
    )

    if "control" in yaml_dict.keys():  # Use a FD control for back exchange correction
        # todo control should be set from an external file
        control_state = yaml_dict["control"]["state"]
        exposure_value = yaml_dict["control"]["exposure"]["value"]
        exposure_units = yaml_dict["control"]["exposure"]["unit"]
        control_exposure = exposure_value * time_factors[exposure_units]

        pmt.set_control((control_state, control_exposure))
    elif (
        "be_percent" in yaml_dict.keys()
    ):  # Flat back exchange percentage for all peptides\
        pmt.set_backexchange(yaml_dict["be_percent"])
    else:
        raise ValueError("No valid back exchange control method specified")

    temperature = yaml_dict["temperature"]["value"]
    try:
        t_offset = temperature_offsets[yaml_dict["temperature"]["unit"]]
    except KeyError:
        t_offset = temperature_offsets[yaml_dict["temperature"]["unit"].lower()]

    temperature += t_offset

    sequence = yaml_dict.get("sequence", "")
    c_term = yaml_dict.get("c_term")
    n_term = yaml_dict.get("n_term") or 1

    if not (c_term or sequence):
        raise ValueError("Must specify either 'c_term' or 'sequence'")

    state_data = pmt.get_state(yaml_dict["state"])
    data_filters = data_filters or []
    for filter in data_filters:
        state_data = filter(state_data)

    hdxm = HDXMeasurement(
        state_data,
        temperature=temperature,
        pH=yaml_dict["pH"],
        sequence=sequence,
        n_term=n_term,
        c_term=c_term,
        **kwargs,
    )

    return hdxm


def load_from_yaml_v040b2(yaml_dict, data_dir=None, **kwargs):  # pragma: no cover
    """
    This is the legacy version to load yaml files of PyHDX v0.4.0b2

    Creates a :class:`~pyhdx.models.HDXMeasurement` object from dictionary input.

    Dictionary can be generated from .yaml format. See templates/yaml_files/SecB.yaml for format specification.

    Parameters
    ----------
    yaml_dict : :obj:`dict`
        Input dictionary specifying experimental metadata and file location to load
    data_dir : :obj:`str` or pathlib.Path object

    Returns
    -------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
        Output data object as specified by `yaml_dict`.
    """

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in yaml_dict["filenames"]]
    else:
        input_files = yaml_dict["filenames"]

    data = read_dynamx(*input_files)

    pmt = PeptideMasterTable(
        data, d_percentage=yaml_dict["d_percentage"]
    )  # todo add proline, n_term options
    if "control" in yaml_dict.keys():  # Use a FD control for back exchange correction
        pmt.set_control(tuple(yaml_dict["control"]))
    elif (
        "be_percent" in yaml_dict.keys()
    ):  # Flat back exchange percentage for all peptides\
        pmt.set_backexchange(yaml_dict["be_percent"])
    else:
        raise ValueError("No valid back exchange control method specified")

    if yaml_dict["temperature_unit"].lower() == "celsius":
        temperature = yaml_dict["temperature"] + 273.15
    elif yaml_dict["temperature_unit"].lower() == "kelvin":
        temperature = yaml_dict["temperature"]
    else:
        raise ValueError(
            "Invalid option for 'temperature_unit', must be 'Celsius' or 'Kelvin'"
        )

    sequence = yaml_dict.get("sequence", "")
    c_term = yaml_dict.get("c_term", 0)
    n_term = yaml_dict.get("n_term", 1)

    if not (c_term or sequence):
        raise ValueError("Must specify either 'c_term' or 'sequence'")

    state_data = pmt.get_state([yaml_dict["series_name"]])
    hdxm = HDXMeasurement(
        state_data,
        temperature=temperature,
        pH=yaml_dict["pH"],
        sequence=sequence,
        n_term=n_term,
        c_term=c_term,
        **kwargs,
    )

    return hdxm
