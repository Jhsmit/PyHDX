from pyhdx.panel.controller import Controller
from pyhdx.support import np_from_txt


def reload_previous(dic, ctrl=None):
    cluster = dic['cluster'] if 'cluster' in dic else None
    ctrl = Controller('template', ['asdf'], cluster=cluster) if ctrl is None else ctrl

    if 'file_paths' not in dic.keys() and 'file_path' in dic.keys():
        dic['file_paths'] = [dic['file_path']]

    while len(ctrl.file_input.file_selectors) < len(dic['file_paths']):
        ctrl._action_add()

    # Read files and add them into the file widgets
    for fs, file_path in zip(ctrl.file_input.file_selectors, dic['file_paths']):
        with open(file_path, 'rb') as f:
            binary = f.read()
        fs.value = binary

    #todo prolines, drop_first

    ctrl.file_input._action_load()

    # Set back exchange parameters
    try:
        ctrl.file_input.norm_mode = dic['norm_mode']
    except KeyError:
        pass

    if ctrl.file_input.norm_mode == 'Exp':
        ctrl.file_input.norm_state = dic['norm_state']
        ctrl.file_input.norm_exposure = dic['norm_exposure']
    else:
        ctrl.file_input.be_percent = dic['be_percent']

    ctrl.file_input.exp_state = dic['exp_state']

    # Apply back exchange correction
    ctrl.file_input._action_parse()

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

    return ctrl

    #todo future make JSON from parameter objects