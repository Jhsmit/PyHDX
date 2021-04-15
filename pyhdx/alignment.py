import re
from pathlib import Path
import numpy as np
import pandas as pd

def parse_clustal_string(s, num_proteins, whitelines=2, offset=0):
    """
    Takes input clustal result and parses it to dictionary with concatenated aligned sequences.

    Parameters
    ----------
    s : :obj:`str`:
        Input clustal string.
    num_proteins : :obj:`int`
        Number of aligned proteins in the clustal result.
    whitelines : :obj:`int`
        Number of white lines between each block of alinged proteins.
    offset : :obj:`int`
        Number of lines before alignment information starts

    Returns
    -------
    alignment : :obj:`dict`
        Concatenated aligned result

    """
    spacing = num_proteins + whitelines
    lines = s.split('\n')
    results = [''.join([re.search('(?<=\s{3})(.*)(?=\t)', line)[0].strip() for line in lines[i + offset::spacing]]) for i in range(num_proteins)]
    names = [re.search('(.*)(?=\s{3})', lines[offset + i])[0].strip() for i in range(num_proteins)]

    alignment = {name: result for name, result in zip(names, results)}

    return alignment


def align_dataframes(dataframes, alignment):
    """
    Aligned dataframes based on an alignment.

    The supplied dataframes should have the residue number as index. The returned dataframes are reindex and their
    residue numbers are moved to the 'r_number' columns.

    Parameters
    ----------
    dataframes : :obj:`list`
        List of input dataframes
    alignment : :obj:`list` or :obj:`dict`
        Alignment as list or dict, values are strings with single-letter amino acids codes where gaps are '-'.

    Returns
    -------

    dataframe: :class:`pd.Dataframe`
        Aligned and concatenated dataframe. If 'alignment' is given as a dict, the returned dataframe is a column
        multiindex dataframe.

    """
    if isinstance(alignment, dict):
        align_list = list(alignment.values())
    else:
        align_list = alignment

    dfs = []
    for df, align in zip(dataframes, align_list):
        align_array = np.array(list(align))
        r_number = np.cumsum(align_array != '-').astype(float)
        r_number[align_array == '-'] = np.nan
        index = pd.Index(r_number, name='r_number')

        result = df.reindex(index).reset_index().astype({'r_number': 'Int32'})
        dfs.append(result)

    keys = alignment.keys() if isinstance(alignment, dict) else None
    output = pd.concat(dfs, axis=1, keys=keys)
    return output