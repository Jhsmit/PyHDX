from __future__ import annotations

import warnings
from typing import Optional, Literal, Union

import pandas as pd
import numpy as np

from pyhdx.support import convert_time, dataframe_intersection


def parse_temperature(value: float, unit: Literal["Celsius", "C", "Kelvin", "K"]):
    temperature_offsets = {"c": 273.15, "celsius": 273.15, "k": 0, "kelvin": 0}

    return value + temperature_offsets[unit.lower()]


COLUMN_ORDER = [
    "start",
    "end",
    "stop",
    "sequence",
    "state",
    "exposure",
    "uptake_corrected",
    "uptake",
    "uptake_sd",
    "maxuptake",
    "ex_residues",
    "fd_uptake",
    "fd_uptake_sd",
    "nd_uptake",
    "nd_uptake_sd",
    "rfu",
    "rfu_sd",
]


def sort_columns(df: pd.DataFrame, column_order: Optional[list[str]] = None) -> pd.DataFrame:
    """
    Sorts columns in DataFrame by a given column order. Columns not in suplied order are appended to the end
    of the dataframe.


    Args:
        df: DataFrame to sort.
        column_order: Order of columns to use. If `None`, a default order is used.

    Returns:
        The column-sorted DataFrame.
    """

    # https://stackoverflow.com/questions/41968732/set-order-of-columns-in-pandas-dataframe
    column_order = column_order or COLUMN_ORDER

    columns_to_order = [col for col in column_order if col in df.columns]
    new_columns = columns_to_order + [col for col in df.columns if col not in columns_to_order]

    return df[new_columns]


def apply_control(
    experiment: pd.DataFrame,
    fd_control: pd.DataFrame,
    nd_control: Optional[pd.DataFrame] = None,
    deepcopy: bool = False,
) -> pd.DataFrame:
    if nd_control is None:
        nd_control = fd_control[["start", "stop", "uptake", "uptake_sd"]].copy()
        nd_control["uptake"] = 0
        nd_control["uptake_sd"] = 0

    intersected = dataframe_intersection(
        [experiment, fd_control, nd_control], ["start", "stop"], reset_index=False
    )

    # select out uptake (u; experiment), FD uptake (f) and ND uptake (n), as well as their sd's
    u, f, n = (
        intersected[0]["uptake"],
        intersected[1]["uptake"],
        intersected[2]["uptake"],
    )
    u_sd, f_sd, n_sd = (
        intersected[0]["uptake_sd"],
        intersected[1]["uptake_sd"],
        intersected[2]["uptake_sd"],
    )

    out_df = intersected[0].copy(deep=deepcopy)
    out_df["rfu"] = (u - n) / (f - n)
    out_df["rfu_sd"] = np.sqrt(
        (1 / (f - n)) ** 2 * u_sd**2
        + ((u - f) / (f - n) ** 2) ** 2 * n_sd**2
        + ((n - u) / (f - n) ** 2) ** 2 * f_sd**2
    )

    out_df["fd_uptake"] = f
    out_df["fd_uptake_sd"] = f_sd
    out_df["nd_uptake"] = n
    out_df["nd_uptake_sd"] = n_sd

    return sort_columns(out_df.reset_index())


def correct_d_uptake(
    peptides: pd.DataFrame,
    drop_first: int = 1,
    d_percentage: float = 100,
    deepcopy: bool = False,
):
    """
    Corrects for back exchange, percentage deuterium in solution and prolines. Adds the number of effective exchanging
    residues as well as corrected deuterium uptake (requires the field 'rfu')

    Modified the 'sequence' field, where n_terminal non-exchanging residues are marked with 'x' and prolines with
    lower case 'p'. Adds the fields '_start' and '_stop', which are the start and stop residue numbers for each peptide
    minus non-exchanging residues.

    Args:
        peptides: DataFrame with peptides
        drop_first: Number of n-terminal residues to consider as fully back-exchanging.
        d_percentage: Percentate deutrium in the exchange buffer.
        deepcopy: Set to `True` to make a deep copy of the input DataFrame, otherwise a shallow copy is made.

    Returns:

    """

    peptides = peptides.copy(deep=deepcopy)

    if not 0.0 <= d_percentage <= 100.0:
        raise ValueError(f"Deuterium percentage must be 0-100, got {d_percentage}")

    peptides["_sequence"] = peptides["sequence"].copy()
    peptides["sequence"] = [s.replace("P", "p") for s in peptides["sequence"]]

    # Find the total number of n terminal / c_terminal residues to remove
    n_term = np.array(
        [len(seq) - len(seq[drop_first:].lstrip("p")) for seq in peptides["sequence"]]
    )

    c_term = np.array([len(seq) - len(seq.rstrip("p")) for seq in peptides["sequence"]])

    peptides["sequence"] = ["x" * nt + s[nt:] for nt, s in zip(n_term, peptides["sequence"])]
    peptides["_start"] = peptides["start"] + n_term
    peptides["_stop"] = peptides["stop"] - c_term

    ex_residues = (
        np.array([len(s) - s.count("x") - s.count("p") for s in peptides["sequence"]])
        * d_percentage
        / 100.0
    )

    peptides["ex_residues"] = ex_residues

    if "rfu" in peptides:
        peptides["uptake_corrected"] = peptides["rfu"] * ex_residues

    return peptides


def verify_sequence(
    df: pd.DataFrame,
    sequence: Optional[str] = None,
    n_term: Optional[int] = None,
    c_term: Optional[int] = None,
) -> tuple[pd.Series, pd.Series]:
    """
    Verify if sequence information in the dataframe is compatible with an externally supplied sequence and/or the residue
    numbers of the N terminal and C terminal residues.

    Args:
        df: Peptide dataframe. Must have columnse 'start', 'stop' and 'sequence'
        sequence: Sequence to check as FASTA string.
        n_term: Optional residue number of N terminal residue. Can be negative to include purification tags.
        c_term: Optional residue number of C terminal residue.

    Returns:
        Tuple of pandas series with full and reconstructed (+ lower case prolines)
    """

    # TODO return single pd series
    # Returns:
    #     Pandas Series with (reconstructed) sequence with residue numbers as index.

    n_term = n_term if n_term is not None else 1

    if sequence is None and c_term is None:
        raise ValueError("Must provide either 'c_term' or 'sequence'")
    elif c_term is None:
        c_term = len(sequence) + n_term - 1

    r_number = pd.RangeIndex(n_term, c_term + 1, name="r_number")

    if df["start"].min() < n_term:
        raise ValueError(
            f"Peptide dataframe contains peptides with start residue number below supplied 'n_term' ({n_term})"
        )
    if df["end"].max() > c_term:
        raise ValueError(
            f"Peptide dataframe contains peptides with end residue number above supplied 'c_term' ({c_term})"
        )

    seq_full = pd.Series(index=r_number, dtype="U").fillna("X")
    seq_reconstruct = pd.Series(index=r_number, dtype="U").fillna("X")

    # iterate over dataframe from C terminal peptides to N terminal peptides
    # paste sequence information in pd.Series at the correct positions.
    for idx in df.index[::-1]:
        start, end = df.loc[idx, "start"], df.loc[idx, "stop"]
        seq_full.loc[start : end - 1] = list(df.loc[idx, "_sequence"])
        seq_reconstruct.loc[start : end - 1] = list(df.loc[idx, "sequence"])

    if sequence:
        for r, s1, s2 in zip(r_number, sequence, seq_full):
            if s2 != "X" and s1 != s2:
                raise ValueError(
                    f"Mismatch in supplied sequence and peptides sequence at residue {r}, expected '{s2}', got '{s1}'"
                )
        if len(sequence) != len(seq_full):
            raise ValueError(
                "Invalid length of supplied sequence. Please check 'n_term' and 'c_term' parameters"
            )
        seq_full = pd.Series(index=r_number, data=list(sequence))

    return seq_full, seq_reconstruct


def filter_peptides(
    df: pd.DataFrame,
    state: Optional[str] = None,
    exposure: Union[dict, float, None] = None,
    query: Optional[list[str]] = None,
    dropna: bool = True,
) -> pd.DataFrame:
    """
    Convenience function to filter a peptides DataFrame.

    Args:
        df: Input :class:`pandas.DataFrame`
        state: Name of protein state to select.
        exposure: Exposure value(s) to select. Exposure is given as a :obj:`dict`, with keys "value" or "values" for
            exposure value, and "unit" for the time unit.
        query: Additional queries to pass to :meth:`pandas.DataFrame.query`.
        dropna: Drop rows with NaN uptake entries.

    Example:
        ::

        d = {"state", "SecB WT apo", "exposure": {"value": 0.167, "unit": "min"}
        filtered_df = filter_peptides(df, **d)

    Returns:

    """

    warnings.warn(
        "`filter_peptides` will be moved to the `hdxms-datasets` package", DeprecationWarning
    )
    if state:
        df = df[df["state"] == state]

    if isinstance(exposure, dict):
        if values := exposure.get("values"):
            values = convert_time(values, exposure.get("unit", "s"), "s")
            df = df[df["exposure"].isin(values)]
        elif value := exposure.get("value"):
            value = convert_time(value, exposure.get("unit", "s"), "s")
            df = df[df["exposure"] == value]
    elif isinstance(exposure, float):
        df = df[df["exposure"] == exposure]

    if query:
        for q in query:
            df = df.query(q)

    if dropna:
        df = df.dropna(subset=["uptake"])

    return df
