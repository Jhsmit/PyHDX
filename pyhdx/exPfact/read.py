import json
import numpy as np
import os


def read_assignments(assignment_file):
    """
    Loads assignments of kints to dexp values.
    :param assignment_file: a file assignments of kints to dexp values.
    :return: 2D numpy array containing assignment information.
    """

    with open(assignment_file, 'r') as f:
        lines = f.readlines()
        # assignments = np.zeros((len(lines), 3))
        assignments =  []
        for ii, line in enumerate(lines):
            assignments.append([int(x) for x in line.strip().split(" ")[:3]])
            # assignments[ii][0] = int(split_line[0])
            # assignments[ii][1] = int(split_line[1])
            # assignments[ii][2] = int(split_line[2])
            # for jj, residue in enumerate(seq):
            #     if jj == 0:
            #         assignments[ii][jj+3] = -1
            #     elif jj == len(seq):
            #         break
            #     elif residue == "P":
            #         assignments[ii][jj+3] = -1
            #         prolines.append(assignments[ii][1]+jj-1)
            #     else:
            #         assignments[ii][jj+3] = calckint(seq[jj - 1], residue, jj-1+assignments[ii][1], temperature, pH)
            #     print("@@@kint",ii+1,jj+int(assignments[ii][1]),assignments[ii][jj+3])

    return np.array(assignments)


def read_seq(seq_file):
    """
    :Loads a peptide sequence form a sequence file
    :param seq_file: file containing the sequence
    :return: list containing the sequence
    """
    with open(seq_file) as f:
        return f.read().strip()


def read_pfact(pfact_file):
    """
    Loads protection factors into a list.
    :param pfact_file: file containing protection factors.
    :return: list
    """
    lines = [line.strip().split() for line in open(pfact_file, 'r').readlines()]
    values = [float(x[1]) for x in lines]
    return np.array(values)


def read_time_points(time_points_file):
    """
    load time points into numpy array
    :param time_points_file:
    :return: numpy array containing time points
    """
    return np.loadtxt(time_points_file)


def read_configuration(config_file):

    config = {}
    with open(config_file, 'r') as f:
        config_json = json.load(f)

        config['base'] =                config_json['base']
        config['assignments'] =         read_assignments(os.path.join(config['base'], config_json['assignments_file']))
        #config['assignments'] = os.path.join(config['base'], config_json['assignments_file'])
        config['dexp'], config['times'] =                read_dexp(os.path.join(config['base'], config_json['dexp_file']))
        config['harmonic_factor'] =     config_json['harmonic_term']
        # config['kint'] =                read_kint(os.path.join(config['base'], config_json['kint_file']))
        config['output'] =              config_json['output_file']
        if config_json['pfact_file']:
            config['pfact'] =           read_pfact(os.path.join(config['base'], config_json['pfact_file']))
        else:
            config['pfact'] = ''
        config['predict'] =             config_json['predict']
        config['random_search'] =       config_json['do_random_search']
        config['random_search_steps'] = config_json['random_search_steps']
        config['time_points'] =         config_json['time_points_file']
        config['tolerance'] =           config_json['tolerance']

        if config_json['seq_file']:
            config['sequence'] = read_seq(os.path.join(config['base'], config_json['seq_file']))
            config['res1'] = 1
            config['resn'] = len(config['sequence'])

        config['weights'] = config_json['weights']
        config['pH'] = config_json['pH']
        config['temperature'] = config_json['temperature']

        return config


def read_dexp(dexp_file):
    """
    Reads dexp values and time points from dexp_file and returns them as numpy arrays
    :param dexp_file:
    :return:
    """
    frag_data = [line.strip().split()[1:] for line in open(dexp_file, 'r').readlines()]
    time_data = [line.strip().split()[:1] for line in open(dexp_file, 'r').readlines()]
    array = []
    for row in frag_data:
        a = [float(x) for x in row]
        array.append(a)
    time_points = np.array([float(x) for x in time_data for x in x])
    return np.array(array).T, time_points
