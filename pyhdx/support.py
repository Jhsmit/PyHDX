import numpy as np
import itertools

def get_reduced_blocks(k_series, max_combine=2, max_join=5):
    #todo arguments of these should be coverage, not series
    block_length = list(k_series.cov.block_length.copy())

    i = 0
    subblock = 0
    while i < len(block_length):
        curr_size = block_length[i]
        if curr_size <= max_combine:
            del block_length[i]
            subblock += curr_size
        elif subblock != 0:
            # if block_length > 1:
            block_length.insert(i, subblock)
            subblock = 0
            i += 1
        else:
            i += 1
    if subblock != 0:
        block_length.append(subblock)

    changes = True
    while changes:
        i = 0
        changes = False
        while i < len(block_length):
            curr_size = block_length[i]
            if curr_size < max_join:
                changes = True
                if i == 0:  # beginning of list
                    block_length[i + 1] += curr_size
                    del block_length[i]
                elif i == len(block_length) - 1:  # end of the list
                    block_length[i - 1] += curr_size
                    del block_length[i]
                else:
                    if block_length[i - 1] < block_length[i + 1]:
                        block_length[i - 1] += curr_size
                    else:
                        block_length[i + 1] += curr_size
                    del block_length[i]
            else:
                i += 1

    return block_length


def get_constant_blocks(k_series, block_size=10, initial_block=5):
    num_repeats = (k_series.cov.prot_len - initial_block) // block_size
    remainder = (k_series.cov.prot_len - initial_block) % block_size

    blocks = [initial_block] + [block_size] * num_repeats
    if remainder:
        blocks += [remainder]

    return blocks


def get_original_blocks(k_series):
    block_length = list(k_series.cov.block_length.copy())
    return block_length


def reduce_inter(args):
    """
    #  https://github.com/brentp/interlap/blob/3c4a5923c97a5d9a11571e0c9ea5bb7ea4e784ee/interlap.py#L224
    # MIT Liscence
    >>> reduce_inter([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]
    >>> reduce_inter([(2, 6), (4, 10)])
    [(2, 10)]
    """
    if len(args) < 2: return args
    args.sort()
    ret = [args[0]]
    for next_i, (s, e) in enumerate(args, start=1):
        if next_i == len(args):
            ret[-1] = ret[-1][0], max(ret[-1][1], e)
            break

        ns, ne = args[next_i]
        if e > ns or ret[-1][1] > ns:
            ret[-1] = ret[-1][0], max(e, ne, ret[-1][1])
        else:
            ret.append((ns, ne))
    return ret


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


def np_from_txt(file_path, delimiter='\t'):
    with open(file_path, 'r') as f:
        header = f.readline()

    if header.startswith('#'):
        names = header[2:].split(delimiter)
    else:
        names = None

    return np.genfromtxt(file_path, dtype=None, names=names, skip_header=1, delimiter=delimiter, encoding=None, autostrip=True)


def try_wrap(coverage, wrap, margin=4):
    """Check for a given coverage if the value of wrap is high enough to not have peptides overlapping within margin"""
    x = np.zeros((wrap, coverage.prot_len + margin))
    wrap_gen = itertools.cycle(range(wrap))
    for i, elem in zip(wrap_gen, coverage.data):
        section = x[i, elem['start']: elem['end'] + 1 + margin]
        if np.any(section):
            return False
        section[:] = 1

    return True


def autowrap(coverage, margin=4):
    """Automatically finds wrap value for coverage to not have overlapping peptides within margin"""
    wrap = 5
    while not try_wrap(coverage, wrap, margin=margin):
        wrap += 5
        if wrap > coverage.prot_len:
            break
    return wrap
