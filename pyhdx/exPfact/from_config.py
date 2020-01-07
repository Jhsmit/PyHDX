from .exPfact import run_config
from .read import read_assignments, read_dexp, read_seq
import os
config = {
    'base': r'',
    'temperature': 300,
    'pH': 7.0,
    'harmonic_factor': 0,
    'random_search_steps': 10000,
    'tolerance': 1e-7,
    'weights': None

}


def make_config(name, dir_in, dir_out, config):
    config['base'] = dir_in

    config['dexp'] = read_dexp(os.path.join(dir_in, name + '.dExp'))
    config['assignments'] = read_assignments(os.path.join(dir_in, name + '.list'))
    config['sequence'] = read_seq(os.path.join(dir_in, name + '.list'))

    config['resl'] = 1
    config['resn'] = len(config['sequence'])

    config['output'] = dir_out

    return config
