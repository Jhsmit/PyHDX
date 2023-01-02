from __future__ import annotations

from pathlib import Path
from typing import Union, Literal, Tuple, List, TextIO

import torch as t
import torch.nn as nn
from io import StringIO, BytesIO
import pandas as pd
from typing.io import IO

import pyhdx
import yaml
import json
import re
import shutil
from datetime import datetime

from importlib import import_module

# Dtype of fields in peptide table data
PEPTIDE_DTYPES = {"start": int, "end": int, "stop": int, "_start": int, "_stop": int}


def read_dynamx(
    filepath_or_buffer: Union[Path[str], str, StringIO],
    time_conversion: Tuple[Literal["h", "min", "s"], Literal["h", "min", "s"]] = (
        "min",
        "s",
    ),
) -> pd.DataFrame:

    """
    Reads DynamX .csv files and returns the resulting peptide table as a pandas DataFrame.

    Args:
        filepath_or_buffer: File path of the .csv file or :class:`~io.StringIO` object.
        time_conversion: How to convert the time unit of the field 'exposure'. Format is ('<from>', <'to'>).
            Unit options are 'h', 'min' or 's'.

    Returns:
        Peptide table as a pandas DataFrame.
    """

    if isinstance(filepath_or_buffer, StringIO):
        hdr = filepath_or_buffer.readline().strip("# \n\t")
        filepath_or_buffer.seek(0)
    else:
        with open(filepath_or_buffer, "r") as f_obj:
            hdr = f_obj.readline().strip("# \n\t")

    names = [name.lower().strip("\r\t\n") for name in hdr.split(",")]
    df = pd.read_csv(filepath_or_buffer, header=0, names=names)

    df.insert(df.columns.get_loc("end") + 1, "stop", df["end"] + 1)

    time_lut = {"h": 3600, "min": 60, "s": 1}
    time_factor = time_lut[time_conversion[0]] / time_lut[time_conversion[1]]

    df["exposure"] *= time_factor
    df.columns = df.columns.str.replace(" ", "_")

    return df


def read_header(file_obj: TextIO, comment: str = "#") -> List[str]:
    header = []

    while True:
        line = file_obj.readline()
        line = line.decode() if isinstance(line, bytes) else line
        if line.startswith(comment):
            header.append(line)
        else:
            break
    return header


def parse_header(filepath_or_buffer, comment="#"):
    if isinstance(filepath_or_buffer, StringIO):
        header = read_header(filepath_or_buffer, comment=comment)
        filepath_or_buffer.seek(0)
    else:
        with open(filepath_or_buffer, "r") as file_obj:
            header = read_header(file_obj, comment=comment)

    header = [h.strip("#\n ") for h in header]
    pattern = r"<[^>]+>"
    header_dict = {}
    for line in header:
        tags = re.findall(r"<[^>]+>", line)
        if len(tags) == 2 and tags[0] == tags[1].replace("/", ""):
            name = tags[0].strip("<>")
            content = json.loads(re.sub(pattern, "", line))
            header_dict[name] = content

    return header_dict


def csv_to_dataframe(filepath_or_buffer, comment="#", **kwargs):
    """
    Reads a .csv file or buffer into a :pandas:`DataFrame` object.
    Comment lines are parsed where json dictionaries marked by tags are read.
    The <pandas_kwargs> marked json dict is used as kwargs for `pd.read_csv`
    The <metadata> marked json dict is stored in the returned dataframe object as `df.attrs['metadata'].

    Parameters
    ----------
    filepath_or_buffer : :obj:`str`, pathlib.Path or io.StringIO
        Filepath or StringIO buffer to read.
    comment : :obj:`str`
        Indicates which lines are comments.
    kwargs
        Optional additional keyword arguments passed to `pd.read_csv`
    Returns
    -------
    df: pd.DataFrame
    """

    if comment is not None:
        header_dict = parse_header(filepath_or_buffer, comment=comment)
    else:
        header_dict = {}

    pd_kwargs = header_dict.get("pandas_kwargs", {})
    pd_kwargs.update(kwargs)
    df = pd.read_csv(filepath_or_buffer, **pd_kwargs)
    if "metadata" in header_dict:
        df.attrs["metadata"] = header_dict["metadata"]
    return df


def csv_to_protein(filepath_or_buffer, comment="#", **kwargs):
    """
    Reads a .csv file or buffer into a  pyhdx.models.Protein object.
    Comment lines are parsed where json dictionaries marked by tags are read.
    The <pandas_kwargs> marked json dict is used as kwargs for `pd.read_csv`
    The <metadata> marked json dict is stored in the returned dataframe object as `df.attrs['metadata'].

    Parameters
    ----------
    filepath_or_buffer : :obj:`str` or pathlib.Path or io.StringIO
        Filepath or StringIO buffer to read.
    comment : :obj:`str`
        Indicates which lines are comments.
    **kwargs : :obj:`dict`, optional
        Optional additional keyword arguments passed to `pd.read_csv`
    Returns
    -------
    protein : pyhdx.models.Protein
        Resulting Protein object with `r_number` as index
    """

    df = csv_to_dataframe(filepath_or_buffer, comment=comment, **kwargs)
    metadata = df.attrs.pop("metadata", {})
    protein = pyhdx.models.Protein(df, **metadata)
    return protein


def csv_to_hdxm(filepath_or_buffer, comment="#", **kwargs):
    """
    Reads a pyhdx .csv file or buffer into a  pyhdx.models.HDXMeasurement or pyhdx.models.HDXMeasurementSet
    object.

    Parameters
    ----------
    filepath_or_buffer : :obj:`str` or pathlib.Path or io.StringIO
        Filepath or StringIO buffer to read.
    comment : :obj:`str`
        Indicates which lines are comments.
    **kwargs : :obj:`dict`, optional
        Optional additional keyword arguments passed to `pd.read_csv`
    Returns
    -------
    protein : pyhdx.models.HDXMeasurement
        Resulting HDXMeasurement object with `r_number` as index
    """

    df = csv_to_dataframe(filepath_or_buffer, comment=comment, **kwargs)
    metadata = df.attrs.pop("metadata", {})
    if df.columns.nlevels == 2:
        hdxm_list = []
        for state in df.columns.unique(level=0):
            subdf = df[state].dropna(how="all").astype(PEPTIDE_DTYPES)
            m = metadata.get(state, {})
            hdxm = pyhdx.models.HDXMeasurement(subdf, **m)
            hdxm_list.append(hdxm)
        data_obj = pyhdx.models.HDXMeasurementSet(hdxm_list)
    elif df.columns.nlevels == 1:
        data_obj = pyhdx.models.HDXMeasurement(df, **metadata)
    else:
        raise ValueError(
            f"Invalid number of column levels, found {df.columns.nlevels}, supported 1 or 2"
        )
    return data_obj


def dataframe_to_stringio(
    df, sio=None, fmt="csv", include_metadata=True, include_version=True, **kwargs
):
    """
    Save a pd.DataFrame to an io.StringIO object. Kwargs to read the resulting .csv object with pd.read_csv to
    get the original pd.DataFrame back are included in the comments.
    Optionally additional metadata or the version of PyHDX used can be included in the comments.

    Parameters
    ----------
    df: pd.DataFrame
        The pandas dataframe to write to the io.StringIO object.
    sio: `io.StringIO`, optional
        The `io.StringIO` object to write to. If `None`, a new `io.StringIO` object is created.
    fmt: :obj:`str`
        Specify the formatting of the output. Options are 'csv' (machine readable) or 'pprint' (human readable)
    include_metadata: :obj:`bool` or :obj:`dict`
        If `True`, the metadata in df.attrs['metadata'] is included. If :obj:`dict`, this dictionary is used as the
        metadata. If `False`, no metadata is included.
    include_version : :obj:`bool`
        `True` to include PyHDX version information.
    **kwargs : :obj:`dict`, optional
            Optional additional keyword arguments passed to `df.to_csv`


    Returns
    -------
    sio: io.StringIO
        Resulting io.StringIO object.

    """
    sio = sio or StringIO()

    if include_version:
        prefix = "# " if fmt == "csv" else ""
        sio.write(prefix + pyhdx.VERSION_STRING + " \n")
        now = datetime.now()
        sio.write(
            prefix + f'{now.strftime("%Y/%m/%d %H:%M:%S")} ({int(now.timestamp())}) \n'
        )

    json_header = {}
    if include_metadata and "metadata" in df.attrs:
        json_header["metadata"] = df.attrs["metadata"]
    elif isinstance(include_metadata, dict):
        json_header["metadata"] = include_metadata

    if fmt == "csv":
        json_header["pandas_kwargs"] = {
            "comment": "#",
            "header": list(range(df.columns.nlevels)),
            "index_col": 0,
        }
        for k, v in json_header.items():
            if v:
                sio.write(f"# <{k}>{json.dumps(v)}</{k}>\n")
        df.to_csv(sio, line_terminator="\n", **kwargs)
    elif fmt == "pprint":
        if include_version:
            sio.write("\n")
        for k, v in json_header.items():
            if v:
                sio.write(f'{k.capitalize().replace("_", " ")}\n')
                sep = len(k) * "-"
                sio.write(f"{sep}\n")
                sio.write(yaml.dump(v, sort_keys=False))
                sio.write("\n")
        # use df.to_string()?
        with pd.option_context(
            "display.max_rows",
            None,
            "display.max_columns",
            None,
            "display.expand_frame_repr",
            False,
        ):
            sio.write(df.__str__())
    else:
        raise ValueError(
            f"Invalid specification for fmt: '{fmt}', must be 'csv' or 'pprint'"
        )

    sio.seek(0)
    return sio


def dataframe_to_file(
    file_path, df, fmt="csv", include_metadata=True, include_version=False, **kwargs
):
    """
    Save a pd.DataFrame to an io.StringIO object. Kwargs to read the resulting .csv object with pd.read_csv to
    get the original pd.DataFrame back are included in the comments.
    Optionally additional metadata or the version of PyHDX used can be included in the comments.

    Parameters
    ----------
    file_path: :obj:`str` or `pathlib.Path`
        File path of the target file to write.
    df: pd.DataFrame
        The pandas dataframe to write to the file.
    fmt: :obj:`str`
        Specify the formatting of the output. Options are '.csv' (machine readable) or 'pprint' (human readable)
    include_metadata: :obj:`bool` or :obj:`dict`
        If `True`, the metadata in df.attrs['metadata'] is included. If :obj:`dict`, this dictionary is used as the
        metadata. If `False`, no metadata is included.
    include_version : :obj:`bool`
        `True` to include PyHDX version information.
    **kwargs : :obj:`dict`, optional
            Optional additional keyword arguments passed to `df.to_csv`


    Returns
    -------
    sio: io.StringIO
        Resulting io.StringIO object.

    """
    sio = dataframe_to_stringio(
        df,
        fmt=fmt,
        include_metadata=include_metadata,
        include_version=include_version,
        **kwargs,
    )
    with open(file_path, "w") as f:
        sio.seek(0)
        shutil.copyfileobj(sio, f)


def save_fitresult(output_dir, fit_result, log_lines=None):
    """
    Save a fit result object to the specified directory with associated metadata

    Output directory contents:
    dG.csv/.txt: Fit output result (dG, covariance, k_obs, pfact)
    losses.csv/.txt: Losses per epoch
    log.txt: Log file with additional metadata (number of epochs, final losses, pyhdx version, time/date)

    Parameters
    ----------
    output_dir: pathlib.Path or :obj:`str`
        Output directory to save fitresult to
    fit_result: pydhx.fittin_torch.TorchFitResult
        fit result object to save
    log_lines: :obj:`list`
        Optional additional lines to write to log file.

    Returns
    -------

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fit_result.to_file(output_dir / "fit_result.csv")
    fit_result.to_file(output_dir / "fit_result.txt", fmt="pprint")

    dataframe_to_file(output_dir / "losses.csv", fit_result.losses)
    dataframe_to_file(output_dir / "losses.txt", fit_result.losses, fmt="pprint")

    if isinstance(
        fit_result.hdxm_set, pyhdx.HDXMeasurement
    ):  # check, but this should always be hdxm_set
        fit_result.hdxm_set.to_file(output_dir / "HDXMeasurement.csv")
    if isinstance(fit_result.hdxm_set, pyhdx.HDXMeasurementSet):
        fit_result.hdxm_set.to_file(output_dir / "HDXMeasurements.csv")

    loss = (
        f"Total_loss {fit_result.total_loss:.2f}, mse_loss {fit_result.mse_loss:.2f}, reg_loss {fit_result.reg_loss:.2f}"
        f"({fit_result.regularization_percentage:.2f}%)"
    )
    epochs = f"Number of epochs: {len(fit_result.losses)}"
    version = pyhdx.VERSION_STRING
    now = datetime.now()
    date = f'# {now.strftime("%Y/%m/%d %H:%M:%S")} ({int(now.timestamp())})'

    lines = [date, version, loss, epochs]
    if log_lines is not None:
        lines.append("")
        lines += log_lines
    log_file_out = output_dir / "log.txt"
    log_file_out.write_text("\n".join(lines))


def load_fitresult(fit_dir):
    """
    Load a fitresult into a fitting_torch.TorchSingleFitResult or :class:`~pyhdx.fitting_torch.TorchBatchFitResult` object

    The fit result must be in the format as generated by saving a fit result with `save_fitresult`.

    Parameters
    ----------
    fir_dir :obj:`str` or :class:`~pathlib.Path`
        Fit result directory.

    """
    pth = Path(fit_dir)
    if pth.is_dir():
        fit_result = csv_to_dataframe(fit_dir / "fit_result.csv")
        losses = csv_to_dataframe(fit_dir / "losses.csv")

        data_obj = csv_to_hdxm(fit_dir / "HDXMeasurements.csv")
        result_klass = pyhdx.fitting_torch.TorchFitResult
    elif pth.is_file():
        raise DeprecationWarning(
            "`load_fitresult` only loads from fit result directories"
        )
        fit_result = csv_to_dataframe(fit_dir)
        assert isinstance(
            hdxm, pyhdx.HDXMeasurement
        ), "No valid HDXMeasurement data object supplied"
    else:
        raise ValueError("Specified fit result path is not a directory")

    fit_metadata = fit_result.attrs.pop("metadata")
    model_klass = getattr(
        import_module("pyhdx.fitting_torch"), fit_metadata["model_name"]
    )

    if isinstance(fit_result.columns, pd.MultiIndex):
        g_arr = fit_result.xs("_dG", level=-1, axis=1).to_numpy().T
    else:
        g_arr = fit_result["_dG"].to_numpy().T
    g_parameter = nn.Parameter(t.tensor(g_arr)).unsqueeze(
        -1
    )  # todo record/generalize shapes
    model = model_klass(g_parameter)

    fit_result_obj = result_klass(data_obj, model, losses=losses, metadata=fit_metadata)

    return fit_result_obj
