import yaml
from pathlib import Path
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.fileIO import read_dynamx
import asyncio
import copy


def load_from_yaml(yaml_dict, data_dir=None, **kwargs):  #name: load what from yaml?
    #todo perhas classmethod on HDXMeasurement object?
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

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in yaml_dict['filenames']]
    else:
        input_files = yaml_dict['filenames']

    data = read_dynamx(*input_files)

    pmt = PeptideMasterTable(data, d_percentage=yaml_dict['d_percentage'])  #todo add proline, n_term options
    if 'control' in yaml_dict.keys():  # Use a FD control for back exchange correction
        pmt.set_control(tuple(yaml_dict['control']))
    elif 'be_percent' in yaml_dict.keys():  # Flat back exchange percentage for all peptides\
        pmt.set_backexchange(yaml_dict['be_percent'])
    else:
        raise ValueError('No valid back exchange control method specified')

    if yaml_dict['temperature_unit'].lower() == 'celsius':
        temperature = yaml_dict['temperature'] + 273.15
    elif yaml_dict['temperature_unit'].lower() == 'kelvin':
        temperature = yaml_dict['temperature']
    else:
        raise ValueError("Invalid option for 'temperature_unit', must be 'Celsius' or 'Kelvin'")

    sequence = yaml_dict.get('sequence', '')
    c_term = yaml_dict.get('c_term', 0)
    n_term = yaml_dict.get('n_term', 1)

    if not (c_term or sequence):
        raise ValueError("Must specify either 'c_term' or 'sequence'")

    state_data = pmt.get_state([yaml_dict['series_name']])
    hdxm = HDXMeasurement(state_data, temperature=temperature, pH=yaml_dict['pH'],
                          sequence=sequence, n_term=n_term, c_term=c_term, **kwargs)

    return hdxm
