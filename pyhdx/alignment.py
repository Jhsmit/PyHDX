from __future__ import annotations

import re
import numpy as np
import pandas as pd


def parse_clustal_string(
    s: str, num_proteins: int, whitelines: int = 2, offset: int = 0
) -> dict:
    """Takes input Clustal result and parses it to a dictionary.

    Keys in the output dict are IDs of the protein as input into clustal. Values are aligned
    (containing '-' for gaps) and concatenated FASTA sequences.

    Args:
        s: Input Clustal string.
        num_proteins: Number of aligned proteins in the clustal result.
        whitelines: Number of white lines between each block of aligned proteins. Default
            value is 2.
        offset: Number of lines before alignment information starts.

    Returns:
        Dictionary with concatenated aligned result

    """
    spacing = num_proteins + whitelines
    lines = s.split("\n")
    results = [
        "".join(
            [
                re.search(r"(?<=\s{3})(.*)(?=\t)", line)[0].strip()
                for line in lines[i + offset :: spacing]
            ]
        )
        for i in range(num_proteins)
    ]
    names = [
        re.search(r"(.*)(?=\s{3})", lines[offset + i])[0].strip()
        for i in range(num_proteins)
    ]

    alignment = {name: result for name, result in zip(names, results)}

    return alignment


# TODO should take dicts only
def align_dataframes(dataframes, alignment, first_r_numbers=None):

    """
    Aligned dataframes based on an alignment.

    The supplied dataframes should have the residue number as index. The returned dataframes are reindex and their
    residue numbers are moved to the 'r_number' columns.

    Parameters
    ----------
    dataframes : :obj:`list` or :obj:`dict`
        Input dataframes
    alignment : :obj:`list` or :obj:`dict`
        Alignment as list or dict, (type must be equal to `dataframes`, values are strings with single-letter amino
        acids codes where gaps are '-'.
    first_r_numbers: :obj:`list`
        List of residue numbers corresponding to the first residue in the alignment sequences. Use for N-terminally
        truncated proteins or proteins with fused purification tags.
    Returns
    -------

    dataframe: :class:`pd.Dataframe`
        Aligned and concatenated dataframe. If 'alignment' is given as a dict, the returned dataframe is a column
        multiindex dataframe.

    """

    # todo add option to merge/include sequence information in output dataframes
    # todo add dict-like input for dataframes (multiindex dataframe)

    assert len(alignment) == len(
        dataframes
    ), "Length of dataframes and alignments does not match"

    if isinstance(alignment, dict):
        if not isinstance(dataframes, dict):
            raise TypeError(
                "'alignment' and 'dataframes' must either both be `dict` or `list`"
            )
        align_list = list(alignment.values())
        df_list = list(dataframes.values())
        assert alignment.keys() == dataframes.keys(), "Keys of input dicts to not match"
    elif isinstance(alignment, list):
        if not isinstance(alignment, list):
            raise TypeError(
                "'alignment' and 'dataframes' must either both be `dict` or `list`"
            )
        align_list = alignment
        df_list = dataframes
    else:
        raise TypeError("Invalid data type for 'alignment'")

    first_r_numbers = (
        first_r_numbers if first_r_numbers is not None else [1] * len(dataframes)
    )
    assert len(first_r_numbers) == len(
        df_list
    ), "Length of first residue number list does not match number of dataframes"

    dfs = []
    for df, align, r_offset in zip(df_list, align_list, first_r_numbers):
        align_array = np.array(list(align))
        r_number = np.cumsum(align_array != "-").astype(float)
        r_number[align_array == "-"] = np.nan
        r_number += r_offset - 1
        index = pd.Index(r_number, name="r_number")

        result = df.reindex(index).reset_index().astype({"r_number": "Int32"})
        dfs.append(result)

    keys = alignment.keys() if isinstance(alignment, dict) else None
    output = pd.concat(dfs, axis=1, keys=keys)
    return output
