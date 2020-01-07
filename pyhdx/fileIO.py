import numpy as np


def write_expfact(data_list, name):
    keyfunc = lambda x: x.exposure
    data_list = sorted(data_list, key=keyfunc)
    check_data(data_list)

    pm = data_list[0]
    write_dexp(data_list, name + '.Dexp')
    write_seq(pm, name + '.seq')
    write_ass(pm, name + '.list')


def check_data(data_list):
    c = data_list[0]
    for d in data_list[1:]:
        for field in ['start', 'end', 'sequence']:
            assert np.all(c.data[field] == d.data[field])


def write_dexp(data_list, fname):
    dexp = np.array([pm.scores / 100 for pm in data_list])
    times = np.array([v.exposure for v in data_list])
    out = np.insert(dexp, 0, times / 60, axis=1)

    np.savetxt(fname, out)


def write_ass(pm, fname):
    new_data = np.empty(len(pm), dtype= [('index', int), ('start', int), ('end', int), ('sequence', '<U64')])
    new_data['index'] = np.arange(len(pm)) + 1
    new_data['start'] = pm.data['start']
    new_data['end'] = pm.data['end']
    new_data['sequence'] = pm.data['sequence']

    #data = np.vstack([np.arange(len(pm)) + 1] + [pm.data[name] for name in ['start', 'end', 'sequence']])
    np.savetxt(fname, new_data, fmt='%i %i %i %10s')


def write_seq(pm, fname):
    with open(fname, 'w') as f:
        f.write(pm.sequence)
