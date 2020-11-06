from pyhdx.panel.main_controllers import PyHDXController
from pyhdx.support import np_from_txt


def reload_previous(dic, ctrl):
    file_input = ctrl.control_panels['PeptideFileInputControl']

    if 'file_paths' not in dic.keys() and 'file_path' in dic.keys():
        dic['file_paths'] = [dic['file_path']]

    while len(file_input.file_selectors) < len(dic['file_paths']):
        ctrl._action_add()

    # Read files and add them into the file widgets
    for fs, file_path in zip(file_input.file_selectors, dic['file_paths']):
        with open(file_path, 'rb') as f:
            binary = f.read()
        fs.value = binary

    #todo prolines, drop_first

    file_input._action_load()

    # Set back exchange parameters
    try:
        file_input.be_mode = dic['be_mode']
    except KeyError:
        pass

    if file_input.be_mode == 'Exp':
        file_input.fd_state = dic['fd_state']
        file_input.fd_exposure = dic['fd_exposure']
    else:
        file_input.be_percent = dic['be_percent']

    file_input.exp_state = dic['exp_state']

    #todo add exp_exposures

    # Apply back exchange correction
    file_input._action_parse()

    if 'sources' not in dic.keys():
        dic['sources'] = {}
    try:
        for k, v in dic['source_files'].items():
            array = np_from_txt(v)
            dic['sources'][k] = {name: array[name] for name in array.dtype.names}
    except KeyError:
        pass

    #todo some of these are generated from fit results
    #there should be a function for that
    for k, v in dic['sources'].items():
        ctrl.publish_data(k, v)

    try:
        for k, v in dic['fit_results'].items():
            ctrl.fit_results[k] = v  #trigger?
    except KeyError:
        pass

    ctrl.param.trigger('fit_results')

    #make this a yaml file
    rcsb_id = dic.get('rcsb_id', '')
    ctrl.control_panels['ProteinViewControl'].rcsb_id = rcsb_id

    return ctrl

    #todo future make JSON from parameter objects
