from __future__ import annotations

import hashlib
from typing import Iterable, List, Mapping, Any

import numpy as np
import itertools
import re
import contextlib
from io import StringIO
from skimage.filters import threshold_multiotsu
import pandas as pd
from itertools import count, groupby
import warnings
from pathlib import Path
from dask.distributed import Client


def make_tuple(item):
    if isinstance(item, list):
        return tuple(make_tuple(i) for i in item)
    elif isinstance(item, dict):
        return tuple((key, make_tuple(value)) for key, value in item.items())
    else:
        return item


def hash_dataframe(df, method="builtin"):
    if method == "builtin":
        tup = (
            *pd.util.hash_pandas_object(df, index=True).values,
            *df.columns,
            *df.columns.names,
            df.index.name,
        )

        return hash(tup)

    elif method == "md5":
        pd_values = list(pd.util.hash_pandas_object(df, index=True).values)
        if isinstance(df.columns, pd.MultiIndex):
            columns = [name for cols in df.columns for name in cols]
        else:
            columns = list(df.columns)

        all_vals = pd_values + columns + list(df.index.names) + list(df.columns.names)
        h = hashlib.md5()
        for val in all_vals:
            h.update(str(val).encode("UTF-8"))

        return h.digest().hex()

    else:
        raise ValueError(f"Invalid method {method!r}, must be 'builtin' or 'md5'")


def hash_array(array, method="builtin"):
    if method == "buitin":
        return hash(array.data.tobytes())
    elif method == "md5":
        h = hashlib.md5(array.data.tobytes())
        return h.digest().hex()
    else:
        raise ValueError(f"Invalid method {method!r}, must be 'builtin' or 'md5'")


def multiindex_apply_function(
    index: pd.MultiIndex,
    level: int,
    func: str,
    args: Iterable = None,
    kwargs: Mapping = None,
) -> pd.MultiIndex:

    args = args or []
    kwargs = kwargs or {}
    new_index = index.set_levels(
        getattr(index.levels[level], func)(*args, **kwargs), level=level
    )

    return new_index


def multiindex_astype(index: pd.MultiIndex, level: int, dtype: str) -> pd.MultiIndex:

    new_index = multiindex_apply_function(index, level, "astype", args=[dtype])
    return new_index


def multiindex_set_categories(
    index: pd.MultiIndex,
    level: int,
    categories: Any,  # index-like
    ordered: bool = False,
    rename: bool = False,
) -> pd.MultiIndex:
    new_index = multiindex_apply_function(
        index,
        level,
        "set_categories",
        args=[categories],
        kwargs=dict(ordered=ordered, rename=rename),
    )
    return new_index


def multiindex_add_categories(
    index: pd.MultiIndex, level: int, categories: Any,  # index-like
) -> pd.MultiIndex:
    new_index = multiindex_apply_function(
        index, level, "add_categories", args=[categories]
    )
    return new_index


def multiindex_astype(index: pd.MultiIndex, level: int, dtype: str,) -> pd.MultiIndex:
    new_index = multiindex_apply_function(index, level, "astype", args=[dtype])
    return new_index


# use tostring(print(df.to_string()))
def df_fullstr(df):
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.expand_frame_repr",
        False,
    ):
        s = df.__str__()

    return s


def pprint_df(df):
    print(df_fullstr(df))


def get_reduced_blocks(coverage, max_combine=2, max_join=5):
    block_length = list(coverage.block_length.copy())

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


def get_constant_blocks(coverage, block_size=10, initial_block=5):
    num_repeats = (coverage.prot_len - initial_block) // block_size
    remainder = (coverage.prot_len - initial_block) % block_size

    blocks = [initial_block] + [block_size] * num_repeats
    if remainder:
        blocks += [remainder]

    return blocks


def get_original_blocks(coverage):
    block_length = list(coverage.block_length.copy())
    return block_length


def reduce_inter(
    args: list[tuple[int, int]], gap_size: int = -1
) -> list[tuple[int, int]]:
    """Reduce overlapping intervals to its non-overlapping intveral parts

    Author: Brent Pedersen
    Source: https://github.com/brentp/interlap/blob/3c4a5923c97a5d9a11571e0c9ea5bb7ea4e784ee/interlap.py#L224

    gap_size : :obj:`int`
        Gaps of this size between adjacent peptides is not considered to overlap. A value of -1 means that peptides
        with exactly zero overlap are separated. With gap_size=0 peptides with exactly zero overlap are not separated,
        and larger values tolerate larger gap sizes.

    >>> reduce_inter([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]
    >>> reduce_inter([(2, 6), (4, 10)])
    [(2, 10)]
    """

    gap_size += 1

    if len(args) < 2:
        return args
    args.sort()
    ret = [args[0]]
    for next_i, (s, e) in enumerate(args, start=1):
        if next_i == len(args):
            ret[-1] = ret[-1][0], max(ret[-1][1], e)
            break

        ns, ne = args[next_i]  # next start, next end
        if (
            e + gap_size > ns or ret[-1][1] + gap_size > ns
        ):  # if current end is further than next start (overlap), OR current inverterval end later then next start
            ret[-1] = (
                ret[-1][0],
                max(e, ne, ret[-1][1]),
            )  # extend the end value of the current inverval by the new end
        else:
            ret.append((ns, ne))
    return ret


@contextlib.contextmanager
def temporary_seed(seed):
    # https://stackoverflow.com/questions/49555991/can-i-create-a-local-numpy-random-seed
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(state)


def grouper(n, iterable, padvalue=None):
    "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
    return itertools.zip_longest(*[iter(iterable)] * n, fillvalue=padvalue)


def _get_f_width(data, sign):
    i = 1 if sign else 0

    w_pos = np.log10(np.nanmax(data)) + i
    w_neg = np.log10(np.nanmax(-data)) + 1

    w = np.nanmax([w_pos, w_neg]) + 1
    try:
        width = int(np.floor(w))
    except OverflowError:  # all zeros
        width = 0
    return width


# move to fileIO?
def fmt_export(
    arr,
    delimiter="\t",
    header=True,
    sig_fig=8,
    width="auto",
    justify="left",
    sign=False,
    pad="",
):
    warnings.warn("fmt_export is to pyhdx.fileIO", PendingDeprecationWarning)
    with np.testing.suppress_warnings() as sup:
        sup.filter(RuntimeWarning)

        flag1 = "" if justify != "left" else "-"
        flag2 = "+" if sign else ""
        flag3 = "0" if pad == "0" else ""
        fmt = []
        hdr = []
        for j, name in enumerate(arr.dtype.names):
            dtype = arr[name].dtype

            if dtype.kind in ["b"]:
                specifier = "i"
                precision = ""
                w = 4 if np.all(arr[name]) else 5
            elif dtype.kind in ["i", "u"]:
                specifier = "i"
                precision = ""

                w = _get_f_width(arr[name], sign)

            elif dtype.kind in ["f", "c"]:
                specifier = "g"
                precision = "." + str(sig_fig)

                # float notation width

                # todo check for nan widths
                w_f = _get_f_width(arr[name], sign) + sig_fig

                # scientific notation width
                i = 1 if sign or np.any(arr[name] < 0) else 0
                w_s = (
                    sig_fig + 4 + i + 1
                )  # +1 for decimal point which is not always needed
                w = min(w_f, w_s) + 1

            elif dtype.kind in ["U", "S", "O"]:
                specifier = "s"
                precision = ""
                w = np.max([len(str(item)) for item in arr[name]])
            else:
                raise TypeError(f"Invalid dtype kind {dtype.kind} for field {name}")

            if width == "auto":
                col_w = w
            elif isinstance(width, int):
                col_w = width
            else:
                raise ValueError("Invalid width")

            if header:
                i = 2 if j == 0 else 0  # Additional space for header comment #
                if width == "auto":
                    _width = max(col_w, len(name) + i)
                elif isinstance(width, int):
                    _width = col_w

                func = str.ljust if justify == "left" else str.rjust
                fill = flag3 if flag3 else " "
                h = func(name, _width - i, fill)
                hdr.append(h)
            else:
                _width = col_w

            s = f"%{flag1}{flag2}{flag3}{_width}{precision}{specifier}"

            fmt.append(s)

    fmt = delimiter.join(fmt)
    hdr = delimiter.join(hdr)
    return fmt, hdr


# move to fileIO?
def np_from_txt(file_path, delimiter="\t"):
    warnings.warn(
        "np_from_txt is moved to pyhdx.fileIO as txt_to_np", PendingDeprecationWarning
    )

    if isinstance(file_path, StringIO):
        file_obj = file_path
    else:
        file_obj = open(file_path, "r")

    names = None
    header_lines = 0
    while True:
        header = file_obj.readline().strip()
        if header.startswith("#"):
            names = header[2:].split(delimiter)
            header_lines += 1
        else:
            break
    file_obj.seek(0)

    return np.genfromtxt(
        file_obj,
        dtype=None,
        names=names,
        skip_header=header_lines,
        delimiter=delimiter,
        encoding=None,
        autostrip=True,
        comments=None,
    )


def try_wrap(start, end, wrap, margin=4):
    """Check for a given coverage if the value of wrap is high enough to not have peptides overlapping within margin

    start, end interval is inclusive, exclusive

    """
    assert len(start) == len(end), "Unequal length of 'start' and 'end' vectors"

    offset = np.min(start)
    start = np.array(start) - offset
    end = np.array(end) - offset

    x = np.zeros((wrap, len(start) + margin))
    wrap_gen = itertools.cycle(range(wrap))
    for i, s, e in zip(wrap_gen, start, end):
        section = x[i, s : e + margin]
        if np.any(section):
            return False
        section[:] = 1

    return True


def autowrap(start, end, margin=4, step=5):
    """
    Automatically finds wrap value for coverage to not have overlapping peptides within margin

    Parameters
    ----------
    start
    end
    margin

    Returns
    -------

    """
    assert len(start) == len(end), "Unequal length of 'start' and 'end' vectors"

    wrap = step
    while True:
        wraps = try_wrap(start, end, wrap, margin=margin)
        wrap += step
        if wraps or wrap > len(start):
            break
    return wrap


# https://stackoverflow.com/questions/15182381/how-to-return-a-view-of-several-columns-in-numpy-structured-array/
def fields_view(arr, fields):
    dtype2 = np.dtype({name: arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)


# https://stackoverflow.com/questions/15182381/how-to-return-a-view-of-several-columns-in-numpy-structured-array/
def make_view(arr, fields, dtype):
    offsets = [arr.dtype.fields[f][1] for f in fields]
    offset = min(offsets)
    stride = max(offsets)
    return np.ndarray(
        (len(arr), 2),
        buffer=arr,
        offset=offset,
        strides=(arr.strides[0], stride - offset),
        dtype=dtype,
    )


vhex = np.vectorize(hex)
base_v = np.vectorize(np.base_repr)


def rgb_to_hex(rgb_a):
    """Converts rgba
    input values are [0, 255]

    alpha is set to zero

    returns as '#000000'

    """
    # Single value
    if isinstance(rgb_a, tuple):
        try:
            r, g, b, a = rgb_a
        except ValueError:
            r, g, b = rgb_a
        return f"#{r:02x}{g:02x}{b:02x}"

    elif isinstance(rgb_a, list):
        try:
            rgba_array = np.array(
                [[b, g, r, 0] for r, g, b, a in rgb_a], dtype=np.uint8
            )
        except ValueError:
            # todo this only works with lists of list and gives to wrong result? tests needed
            rgba_array = np.array([[b, g, r, 0] for r, g, b in rgb_a], dtype=np.uint8)

    elif isinstance(rgb_a, np.ndarray):
        # todo: allow rgb arrays
        assert rgb_a.shape[-1] == 4
        if rgb_a.data.c_contiguous:
            # todo check for c-contigious
            rgba_array = rgb_a
        else:
            rgba_array = np.array(rgb_a)
    else:
        raise TypeError(f"Invalid type for 'rgb_a': {rgb_a}")

    ints = rgba_array.astype(np.uint8).view(dtype=np.uint32).byteswap()
    padded = np.char.rjust(base_v(ints // 2 ** 8, 16), 6, "0")
    result = np.char.add("#", padded).squeeze()

    return result


# def rgb_to_hex(r, g, b):
#     return f'#{r:02x}{g:02x}{b:02x}'


def hex_to_rgb(h):
    """returns rgb as int 0-255"""
    r, g, b = tuple(int(h.lstrip("#")[2 * i : 2 * i + 2], 16) for i in range(3))
    return r, g, b


def hex_to_rgba(h, alpha=255):
    r, g, b = tuple(int(h.lstrip("#")[2 * i : 2 * i + 2], 16) for i in range(3))
    return r, g, b, alpha


def group_with_index(arr):
    # https://stackoverflow.com/questions/25438491/finding-consecutively-repeating-strings-in-python-list/25438531#25438531
    i = 0
    for k, vs in itertools.groupby(arr):
        c = sum(1 for _ in vs)
        yield k, c, i
        i += c


# move to output?
# pymol coloring functions need cleaning/refactoring
# remvoe color_to_pymol
# instead use apply_cmap then color_pymol / color_pymol_script


def colors_to_pymol(r_number, color_arr, c_term=None, no_coverage="#8c8c8c"):
    """coverts colors (hexadecimal format) and corresponding residue numbers to pml
    script to color structures in pymol
    residue ranges in output are inclusive, incluive

    c_term:
        optional residue number of the c terminal of the last peptide doedsnt cover the c terminal
    """

    # todo replace with pandas dataframe magic

    c_term = c_term or np.max(r_number)
    pd_series = pd.Series(color_arr, index=r_number)
    pd_series = pd_series.reindex(np.arange(1, c_term + 1))
    pd_series = pd_series.replace("nan", no_coverage)  # No coverage at nan entries
    pd_series = pd_series.replace(np.nan, no_coverage)  # Numpy NaNs

    return series_to_pymol(pd_series)


def apply_cmap(pd_series_or_df, cmap, norm=None):
    values = pd_series_or_df if norm is None else norm(pd_series_or_df)
    rgb_colors = cmap(values, bytes=True)
    hex_colors = rgb_to_hex(rgb_colors)

    if isinstance(pd_series_or_df, pd.Series):
        return pd.Series(hex_colors, index=pd_series_or_df.index)
    elif isinstance(pd_series_or_df, pd.DataFrame):
        return pd.DataFrame(
            hex_colors, index=pd_series_or_df.index, columns=pd_series_or_df.columns
        )


def color_pymol(pd_series, cmd, model=None):
    grp = pd_series.groupby(pd_series)

    for c, pd_series in grp:
        result = [
            list(g)
            for _, g in groupby(pd_series.index, key=lambda n, c=count(): n - next(c))
        ]
        r, g, b = hex_to_rgb(c)
        residues = [f"resi {g[0]}-{g[-1]}" for g in result]
        selection = " + ".join(residues)

        if model:
            selection = f"model {model} and ({selection})"

        cmd.set_color(c, [r, g, b])
        cmd.color(c, selection=selection)


def series_to_pymol(pd_series):
    """
    Coverts a pandas series to pymol script to color proteins structures in pymol
    Series must have hexadecimal color values and residue number as index

    Parameters
    ----------
    pd_series : :class:`~pandas.Series`

    Returns
    -------

    s_out : :obj:`str`

    """

    # https://stackoverflow.com/questions/33483670/how-to-group-a-series-by-values-in-pandas
    grp = pd_series.groupby(pd_series)

    s_out = ""
    for c, pd_series in grp:
        r, g, b = hex_to_rgb(c)
        s_out += f"set_color color_{c}, [{r},{g},{b}]\n"

    # https://stackoverflow.com/questions/30993182/how-to-get-the-index-range-of-a-list-that-the-values-satisfy-some-criterion-in-p
    for c, pd_series in grp:
        result = [
            list(g)
            for _, g in groupby(pd_series.index, key=lambda n, c=count(): n - next(c))
        ]
        residues = [f"resi {g[0]}-{g[-1]}" for g in result]

        s_out += f"color color_{c}, " + " + ".join(residues) + "\n"

    return s_out


def make_monomer(input_file, output_file):
    """reads input_file pdb file and removes all chains except chain A and all water"""
    with open(input_file, "r") as f_in:
        with open(output_file, "w") as f_out:
            for line in iter(f_in.readline, ""):
                if line.startswith("COMPND") and "CHAIN" in line:
                    res = re.findall(":(.*);", line)[0]
                    line = line.replace(res + ";", " A;" + " " * (len(res) - 2))
                if line.startswith("ATOM") and not " A " in line:
                    continue
                elif line.startswith("HETATM") and "HOH" in line:
                    continue
                f_out.write(line)


# move t
# o output?
def make_color_array(rates, colors, thds, no_coverage="#8c8c8c"):
    """

    :param rates: array of rates
    :param colors: list of colors (slow to fast)
    :param thds: list of thresholds
    no_coverage: color value for no coverage
    :return:
    """

    output = np.full_like(rates, fill_value=no_coverage, dtype="U7")
    full_thds = [-np.inf] + list(thds) + [np.inf]
    for lower, upper, color in zip(full_thds[:-1], full_thds[1:], colors):
        b = (rates > lower) & (rates <= upper)

        output[b] = color

    return output


def multi_otsu(*rates, classes=3):
    """
    global otsu thesholding of multiple rate arrays in log space

    Parameters
    ----------
    rates : iterable
        iterable of numpy structured arrays with  a 'rate' field
    classes : :obj:`int`
        Number of classes to divide the data into

    Returns
    -------
    thds : :obj:`tuple`
        tuple with thresholds

    """
    all_rates = np.concatenate([data["rate"] for data in rates])
    thd_rates = np.log(all_rates[~np.isnan(all_rates)])
    thds = threshold_multiotsu(thd_rates, classes=classes)
    return tuple(np.e ** thd for thd in thds)


def scale(x, out_range=(-1, 1)):
    """rescale input array x to range `out_range`"""
    domain = np.nanmin(x), np.nanmax(x)
    y = (x - (domain[1] + domain[0]) / 2) / (domain[1] - domain[0])
    return y * (out_range[1] - out_range[0]) + (out_range[1] + out_range[0]) / 2


def gen_subclasses(cls):
    """Recursively find all subclasses of cls"""
    for sub_cls in cls.__subclasses__():
        yield sub_cls
        yield from gen_subclasses(sub_cls)


def pprint_df_to_file(df, file_path_or_obj):
    """
    Pretty print (human-readable) a dataframe to a file

    Parameters
    ----------
    df : :class:`~pandas.DataFrame`
    file_path_or_obj : :obj:`str`, Path or :class:`~io.StringIO`

    """
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.expand_frame_repr",
        False,
    ):  # more options can be specified also
        if isinstance(file_path_or_obj, str):
            pth = Path(file_path_or_obj)
            pth.write_text(df.__str__())
        elif isinstance(file_path_or_obj, Path):
            file_path_or_obj.write_text(df.__str__())
        elif isinstance(file_path_or_obj, StringIO):
            file_path_or_obj.write(df.__str__())
