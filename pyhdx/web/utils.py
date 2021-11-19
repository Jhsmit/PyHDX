from pathlib import Path

#todo merge with batch_processing
time_factors = {"s": 1, "m": 60., "min": 60., "h": 3600, "d": 86400}
temperature_offsets = {'c': 273.15, 'celsius': 273.15, 'k': 0, 'kelvin': 0}

def load_state(ctrl, yaml_dict, data_dir, name=None):
    """
    Load a HDXMeasurement into the web interface
    Experimental use only

    Parameters
    ----------
    ctrl: :class:`~pyhdx.web.main_controllers.PyHDXController`
        Main controller to load the data into
    yaml_dict: :obj:`dict`
        Dictionary with measurement info according to batch processing format
    data_dir: path_like
    name: :obj:`str`, optional
        Optional name for the HDXMeasurement

    Returns
    -------
        None
    """
    #raise DeprecationWarning("Currently not up-to-date")

    if data_dir is not None:
        input_files = [Path(data_dir) / fname for fname in yaml_dict['filenames']]
    else:
        input_files = [Path(p) for p in yaml_dict['filenames']]

    files = [f.read_bytes() for f in input_files]

    file_input = ctrl.control_panels['PeptideFileInputControl']
    file_input.input_files = files

    control_state = yaml_dict['control']['state']
    control_exp = yaml_dict['control']['exposure']['value']*time_factors[yaml_dict['control']['exposure']['unit']]

    file_input.fd_state = control_state
    file_input.fd_exposure = control_exp
    file_input.pH = yaml_dict['pH']
    file_input.temperature = yaml_dict['temperature']['value'] + temperature_offsets[yaml_dict['temperature']['unit'].lower()]
    file_input.d_percentage = yaml_dict['d_percentage']

    file_input.exp_state = yaml_dict['state']
    file_input.dataset_name = name or yaml_dict['state']
    file_input._action_add_dataset()