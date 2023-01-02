from typing import Any

from pyhdx import HDXMeasurement, PeptideMasterTable
from pyhdx.batch_processing import StateParser, time_factors, temperature_offsets


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


LOADER_VERSIONS = {"041": load_hdxm_v041}
