import yaml
from pathlib import Path
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.fileIO import read_dynamx, txt_to_protein, csv_to_protein
import asyncio
import copy


def load_from_yaml(yaml_dict, data_dir=None):  #name: load what from yaml?
    #todo perhas classmethod on HDXMeasurement object?
    """
    Creates a :class:`~pyhdx.models.HDXMeasurement` object from dictionary input.

    Dictionary can be generated from .yaml format and should specify

    Parameters
    ----------
    yaml_dict : :obj:`dict`
        Input dictionary specifying metadata and file location to load
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

    try:
        c_term = yaml_dict.get('c_term', 0) or len(yaml_dict['sequence']) + 1
    except KeyError:
        raise ValueError("Must specify either 'c_term' or 'sequence'")
    #todo parse n_term and sequence

    if yaml_dict['temperature_unit'].lower() == 'celsius':
        temperature = yaml_dict['temperature'] + 273.15
    elif yaml_dict['temperature_unit'].lower() == 'kelvin':
        temperature = yaml_dict['temperature']
    else:
        raise ValueError("Invalid option for 'temperature_unit', must be 'Celsius' or 'Kelvin'")

    sequence = yaml_dict.get('sequence', None)

    state_data = pmt.get_state([yaml_dict['series_name']])
    hdxm = HDXMeasurement(state_data, c_term=c_term, temperature=temperature, pH=yaml_dict['pH'], sequence=sequence)

    return hdxm


def load_folding_from_yaml(yaml_dict, data_dir=None):  #name: load what from yaml?
    """


    """

    raise NotImplementedError('Loading folding data from yaml currently not implemented')

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in yaml_dict['filenames']]
    else:
        input_files = yaml_dict['filenames']

    data = read_dynamx(*input_files)

    pmt = PeptideMasterTable(data, d_percentage=yaml_dict['d_percentage'])  #todo add proline, n_term options
    #todo merge this with the other func where it checks for control names to determine what to apply
    pmt.set_control(control_1=tuple(yaml_dict['control_1']), control_0=tuple(yaml_dict['control_0']))

    try:
        c_term = yaml_dict.get('c_term', 0) or len(yaml_dict['sequence']) + 1
    except KeyError:
        raise ValueError("Must specify either 'c_term' or 'sequence'")

    states = pmt.groupby_state(c_term=c_term)
    series = states[yaml_dict['series_name']]

    if yaml_dict['temperature_unit'].lower() == 'celsius':
        temperature = yaml_dict['temperature'] + 273.15
    elif yaml_dict['temperature_unit'].lower() == 'kelvin':
        temperature = yaml_dict['temperature']
    else:
        raise ValueError("Invalid option for 'temperature_unit', must be 'Celsius' or 'Kelvin'")

    kf = KineticsFitting(series, temperature=temperature, pH=yaml_dict['pH'])

    return kf


def do_fitting_from_yaml(yaml_dict, kf_obj):
    raise NotImplementedError('Fitting from yaml not implemented')
    yaml_dict = copy.deepcopy(yaml_dict)  # Make local copy to not affect the supplied dict by pop
    guess = yaml_dict['initial_guess']
    if 'file_path' in guess.keys():
        try:  # todo update when txt files are no longer in existence
            initial_guess = txt_to_protein(guess['file_path'])
        except KeyError:
            initial_guess = csv_to_protein(guess['file_path'])
    else:
        raise NotImplementedError('only guesses by file currently')

    global_fit = yaml_dict['global_fit']
    optimizer_kwargs = global_fit.pop('optimizer_kwargs')
    fit_result = kf_obj.global_fit(initial_guess, **global_fit, **optimizer_kwargs)

    return fit_result
