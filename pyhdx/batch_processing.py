import warnings
from pathlib import Path
import os
from pyhdx.models import PeptideMasterTable, HDXMeasurement, HDXMeasurementSet
from pyhdx.fileIO import read_dynamx


time_factors = {"s": 1, "m": 60.0, "min": 60.0, "h": 3600, "d": 86400}
temperature_offsets = {"c": 273.15, "celsius": 273.15, "k": 0, "kelvin": 0}

# todo add data filters in yaml spec
# todo add proline, n_term options
class YamlParser(object):
    ""'object used to parse yaml data input files into PyHDX HDX Measurement object'

    def __init__(self, yaml_dict, data_src=None, data_dict=None):
        self.yaml_dict = yaml_dict
        if isinstance(data_src, (os.PathLike, str)):
            self.data_src = Path(data_src)
        elif isinstance(data_src, dict):
            self.data_src = data_src
        else:
            raise TypeError(f"Invalid data type {type(data_src)!r}, must be path or dict")

    def load_data(self, *filenames, reader='dynamx'):
        if reader == 'dynamx':
            read_func = read_dynamx
        else:
            raise NotImplementedError("Only reading of dynamx files is implemented")

        if isinstance(self.data_src, Path):
            input_files = [self.data_src / filename for filename in filenames]
            df = read_func(*input_files)
        else:
            input_files = [self.data_src[filename] for filename in filenames]
            df = read_func(*input_files)

        return df

    def load_hdxmset(self):
        """batch read the full yaml spec into a hdxmeasurementset"""
        hdxm_list = []
        for state in self.yaml_dict.keys():
            hdxm = self.load_hdxm(state, name=state)
            hdxm_list.append(hdxm)

        return HDXMeasurementSet(hdxm_list)

    def load_hdxm(self, state, **kwargs):
        """read a single protein state to hdxmeasurement
        kwargs: additional kwargs passed to hdxmeasurementset
        """

        state_dict = self.yaml_dict[state]

        filenames = state_dict["filenames"]
        df = self.load_data(*filenames)

        pmt = PeptideMasterTable(df,
                                 drop_first=state_dict.get('drop_first', 1),
                                 d_percentage=state_dict['d_percentage'])

        if 'control' in state_dict.keys():  # Use a FD control for back exchange correction
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

        hdxm = HDXMeasurement(
            state_data,
            temperature=temperature,
            pH=state_dict["pH"],
            sequence=sequence,
            n_term=n_term,
            c_term=c_term,
            **kwargs
        )

        return hdxm


def yaml_to_hdxmset(yaml_dict, data_dir=None, **kwargs):
    """reads files according to `yaml_dict` spec from `data_dir into HDXMEasurementSet"""

    hdxm_list = []
    for k, v in yaml_dict.items():
        hdxm = yaml_to_hdxm(v, data_dir=data_dir, name=k)
        hdxm_list.append(hdxm)

    return HDXMeasurementSet(hdxm_list)


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

    warnings.warn('This method is deprecated in favor of YamlParser', DeprecationWarning)

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in yaml_dict["filenames"]]
    else:
        input_files = yaml_dict["filenames"]

    data = read_dynamx(*input_files)

    pmt = PeptideMasterTable(data,
                             drop_first=yaml_dict.get('drop_first', 1),
                             d_percentage=yaml_dict['d_percentage'])

    if 'control' in yaml_dict.keys():  # Use a FD control for back exchange correction
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
        **kwargs
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
        **kwargs
    )

    return hdxm
