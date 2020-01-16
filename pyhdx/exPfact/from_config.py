from .read import read_assignments, read_dexp, read_seq
import numpy as np
import os
CONFIG = {
    'base': r'',
    'temperature': 300,
    'pH': 7.0,
    'harmonic_factor': 0,
    'random_search_steps': 10000,
    'tolerance': 1e-7,
    'weights': None,
    'pfact': ''
}


def make_config(name, dir_in, dir_out):
    config = CONFIG.copy()

    config['base'] = dir_in

    config['dexp'], config['times'] = read_dexp(os.path.join(dir_in, name + '.dExp'))

    #config['dexp'] = read_dexp(os.path.join(dir_in, name + '.dExp'))
    config['assignments'] = read_assignments(os.path.join(dir_in, name + '.list'))
    config['sequence'] = read_seq(os.path.join(dir_in, name + '.seq'))

    # full_dex = np.genfromtxt(os.path.join(dir_in,  name + '.dExp'))
    # times = full_dex[:, 0]
    # config['times'] = times
    config['res1'] = 1
    config['resn'] = len(config['sequence'])

    config['output'] = os.path.join(dir_out, name + 'output')

    return config
