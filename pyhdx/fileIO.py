import numpy as np
from numpy.lib.recfunctions import stack_arrays
from io import StringIO
import pandas as pd
import pyhdx

def read_dynamx(*file_paths, intervals=('inclusive', 'inclusive'), time_unit='min'):
    """
    Reads a dynamX .csv file and returns the data as a numpy structured array

    Parameters
    ----------
    file_paths : :obj:`iterable`
        File path of the .csv file or StringIO object
    intervals : :obj:`tuple`
        Format of how start and end intervals are specified.
    time_unit : :obj:`str`
        Not implemented

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
    if intervals[0] == 'inclusive':
        start_correction = 0
    elif intervals[0] == 'exclusive':
        start_correction = 1
    else:
        raise ValueError(f"Invalid start interval value {intervals[0]}, must be 'inclusive' or 'exclusive'")
    if intervals[1] == 'inclusive':
        end_correction = 1
    elif intervals[1] == 'exclusive':
        end_correction = 0
    else:
        raise ValueError(f"Invalid start interval value {intervals[1]}, must be 'inclusive' or 'exclusive'")

    full_data['start'] += start_correction
    full_data['end'] += end_correction

    return full_data


def csv_to_np(file_path, delimiter='\t', column_depth=None):
    """Read csv file and returns a numpy ndarray"""
    if isinstance(file_path, StringIO):
        file_obj = file_path
    else:
        file_obj = open(file_path, 'r')

    lines = []
    num_comments_lines = 0
    while True:
        line = file_obj.readline()
        num_comments_lines += 1
        lines.append(line)
        if not line.startswith('#'):
            break
          #first two lines

    if column_depth is None:
        column_depth = 0

        while True:
            line = file_obj.readline()
            column_depth += 1

            if np.mean([c.isdigit() for c in line]) > 0.1:
                break

            lines.append(line)
    else:
        for i in range(column_depth - 1):
            line = file_obj.readline()
            lines.append(line)

    names = lines[-1].split(',')

    if column_depth > 1:
        raise ValueError('Cannot read multi-index column files into numpy structured array')

    file_obj.seek(0)

    return np.genfromtxt(file_obj, dtype=None, names=names, skip_header=num_comments_lines+column_depth, delimiter=delimiter,
                         encoding=None, autostrip=True, comments=None, deletechars='')


def txt_to_np(file_path, delimiter='\t'):
    """Read .txt file and returns a numpy ndarray"""
    if isinstance(file_path, StringIO):
        file_obj = file_path
    else:
        file_obj = open(file_path, 'r')

    names = None
    header_lines = 0
    while True:
        header = file_obj.readline().strip()
        if header.startswith('#'):
            names = header[2:].split(delimiter)
            header_lines += 1
        else:
            break
    file_obj.seek(0)

    return np.genfromtxt(file_obj, dtype=None, names=names, skip_header=header_lines, delimiter=delimiter,
                         encoding=None, autostrip=True, comments=None, deletechars='')


def txt_to_pd(file_path):

    pass

def parse_header(file_path, comment='#'):
    """
    Parse a the header of a txt file and determine the number of comment and header lines.

    Parameters
    ----------
    file_path

    Returns
    -------

    """
    if isinstance(file_path, StringIO):
        file_obj = file_path
    else:
        file_obj = open(file_path, 'r')

    num_comments_lines = 0
    while file_obj.readline().startswith(comment):
        num_comments_lines += 1

    column_depth = 1
    while np.mean([c.isdigit() for c in file_obj.readline()]) < 0.1:
        column_depth += 1

    file_obj.seek(0)
    return num_comments_lines, column_depth


def csv_to_dataframe(file_path, column_depth=None, **kwargs):
    #todo @tejas: intersphinx + update docstring
    """
    Read .csv file and return <pandas dataframe>

    Parameters
    ----------
    file_path
    column_depth: number of lines after header which describe column headings (typical = 1)

    Returns
    -------

    """
    if isinstance(file_path, StringIO):
        file_obj = file_path
    else:
        file_obj = open(file_path, 'r')

    num_comments_lines = 0
    while file_obj.readline().startswith('#'):
        num_comments_lines += 1

    if column_depth is None:
        column_depth = 1
        while np.mean([c.isdigit() for c in file_obj.readline()]) < 0.1:
            column_depth += 1

    file_obj.seek(0)
    header = [i + num_comments_lines for i in range(column_depth)]
    df = pd.read_csv(file_path, header=header, **kwargs)
    return df


def csv_to_protein(file_path, column_depth=None):
    #todo @tejas: intersphinx + update docstring
    """

    Parameters
    ----------
    file_path
    column_depth

    Returns
    -------

    """

    df = csv_to_dataframe(file_path, column_depth=column_depth, index_col=0)

    return pyhdx.models.Protein(df) #todo metadata


def txt_to_protein(file_path):
    """read .txt file and returns a :class:`pyhdx.models.Protein` object"""
    array = txt_to_np(file_path)
    return pyhdx.models.Protein(array, index='r_number')


def fmt_export(arr, delimiter='\t', header=True, sig_fig=8, width='auto', justify='left', sign=False, pad=''):
    with np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning)

        flag1 = '' if justify != 'left' else '-'
        flag2 = '+' if sign else ''
        flag3 = '0' if pad == '0' else ''
        fmt = []
        hdr = []
        for j, name in enumerate(arr.dtype.names):
            dtype = arr[name].dtype

            if dtype.kind in ['b']:
                specifier = 'i'
                precision = ''
                w = 4 if np.all(arr[name]) else 5
            elif dtype.kind in ['i', 'u']:
                specifier = 'i'
                precision = ''

                w = _get_f_width(arr[name], sign)

            elif dtype.kind in ['f', 'c']:
                specifier = 'g'
                precision = '.' + str(sig_fig)

                # float notation width

                # todo check for nan widths
                w_f = _get_f_width(arr[name], sign) + sig_fig

                # scientific notation width
                i = 1 if sign or np.any(arr[name] < 0) else 0
                w_s = sig_fig + 4 + i + 1  # +1 for decimal point which is not always needed
                w = min(w_f, w_s) + 1

            elif dtype.kind in ['U', 'S', 'O']:
                specifier = 's'
                precision = ''
                w = np.max([len(str(item)) for item in arr[name]])
            else:
                raise TypeError(f'Invalid dtype kind {dtype.kind} for field {name}')

            if width == 'auto':
                col_w = w
            elif isinstance(width, int):
                col_w = width
            else:
                raise ValueError('Invalid width')

            if header:
                i = 2 if j == 0 else 0  # Additional space for header comment #
                if width == 'auto':
                    _width = max(col_w, len(name) + i)
                elif isinstance(width, int):
                    _width = col_w

                func = str.ljust if justify == 'left' else str.rjust
                fill = flag3 if flag3 else ' '
                h = func(name, _width - i, fill)
                hdr.append(h)
            else:
                _width = col_w

            s = f'%{flag1}{flag2}{flag3}{_width}{precision}{specifier}'

            fmt.append(s)

    fmt = delimiter.join(fmt)
    hdr = delimiter.join(hdr)
    return fmt, hdr


def _get_f_width(data, sign):
    i = 1 if sign else 0


    w_pos = np.log10(np.nanmax(data)) + i
    w_neg = np.log10(np.nanmax(-data)) + 1

    w = np.nanmax([w_pos, w_neg]) + 1
    try:
        width = int(np.floor(w))
    except OverflowError: # all zeros
        width = 0
    return width
