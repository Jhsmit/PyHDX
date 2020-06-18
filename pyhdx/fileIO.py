import numpy as np
from numpy.lib.recfunctions import stack_arrays
from io import StringIO


def read_dynamx(*file_paths):
    """
    Reads a dynamX .csv file and returns the data as a numpy structured array

    Parameters
    ----------
    file_paths: :obj:`iterable`
        File path of the .csv file or StringIO object

    Returns
    -------
    data : :class:`~numpy.ndarray`
        Peptides as a numpy structured array

    """

    data_list = []
    for fpath in file_paths:
        # names = [t[0] for t in CSV_DTYPE]
        if isinstance(fpath, StringIO):
            hdr = fpath.readline().strip('# \n\t')
            fpath.seek(0)
        else:
            with open(fpath, 'r') as f:
                hdr = f.readline().strip('# \n\t')

        names = [name.lower() for name in hdr.split(',')]
        data = np.genfromtxt(fpath, skip_header=1, delimiter=',', dtype=None, names=names, encoding='UTF-8')
        data_list.append(data)

    full_data = stack_arrays(data_list, usemask=True, autoconvert=True)
    return full_data

