import yaml
from pathlib import Path
from pyhdx.models import PeptideMasterTable
from pyhdx.fitting import KineticsFitting
from pyhdx.fileIO import read_dynamx, txt_to_protein



def load_from_yaml(yaml_dict, data_dir=None):  #name: load what from yaml?
    """
    Creates a :class:`~pyhdx.fitting.KineticsFitting` object from dictionary input.

    Dictionary can be generated from .yaml format and should specifiy

    Parameters
    ----------
    yaml_dict : :obj:`dict`
        Input dictionary specifying metadata and file location to load
    data_dir : :obj:`str` or pathlib.Path object

    Returns
    -------
    kf : :class:`~pyhdx.fititng.KineticsFitting`
        :class:`~pyhdx.fititng.KineticsFitting` class as specified by input dict.
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
    guess = yaml_dict['initial_guess']
    if 'file_path' in guess.keys():
        initial_guess = txt_to_protein(guess['file_path'])
    else:
        raise NotImplementedError('only guesses by file currently')

    global_fit = yaml_dict['global_fit']
    optimizer_kwargs = global_fit.pop('optimizer_kwargs')
    fit_result = kf_obj.global_fit_torch(initial_guess, **global_fit, **optimizer_kwargs)

    return fit_result


def do_stuff(load_dict, fit_dict, output_dir):
    kf = load_from_yaml(load_dict)
    fr = do_fitting_from_yaml(fit_dict, kf)
    output_dir = Path(output_dir)
    combined_dict = {'load': load_dict, 'fitting': fit_dict}

    yaml_file_out = output_dir / 'settings.yaml'
    yaml_file_out.write_text(yaml.dump(combined_dict))

    fit_file_out = output_dir / 'deltaG.txt'
    fr.output.to_file(fit_file_out)

    loss_file_out = output_dir / 'reg_loss.txt'
    np.savetxt(loss_file_out, fr.metadata['reg_loss'])

    pickle_file_out = output_dir / 'fit_result.pickle'
    with pickle_file_out.open(mode='wb') as f:
        pickle.dump(fr, f)
