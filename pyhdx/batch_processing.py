from __future__ import annotations

import warnings
from dataclasses import dataclass
from io import StringIO
from functools import reduce, cached_property
from pathlib import Path
import os
import re
from typing import Union, Any, Optional, Iterable, Literal

from pyhdx import TorchFitResult
from pyhdx.models import PeptideMasterTable, HDXMeasurement, HDXMeasurementSet
from pyhdx.fileIO import read_dynamx, csv_to_dataframe, save_fitresult
from pyhdx.process import correct_d_uptake, apply_control
from pyhdx.support import gen_subclasses, filter_peptides
from pyhdx.config import cfg

import param
import pandas as pd
import numpy as np
import yaml

time_factors = {"s": 1, "m": 60.0, "min": 60.0, "h": 3600, "d": 86400}
temperature_offsets = {"c": 273.15, "celsius": 273.15, "k": 0, "kelvin": 0}

@dataclass(frozen=True)
class DataFile(object):

    name: str

    format: Literal['DynamX']

    filepath_or_buffer: Union[Path, StringIO]

    @cached_property
    def data(self):
        if self.format == 'DynamX':
            return read_dynamx(self.filepath_or_buffer)

        if isinstance(self.filepath_or_buffer, StringIO):
            self.filepath_or_buffer.seek(0)


class StateParser(object):
    """

    Args:
        hdx_spec: Dictionary with HDX-MS state specification.
        data_src: Optional data source with input data files. If not specified, current
            directory is used. Otherwise, either a data source path can be specified or
            data can be given as a dictionary, where keys are filenames and values are
            :class:`~io.StringIO` with file contents.
    """

    def __init__(
        self,
        hdx_spec: dict,
        data_src: Union[os.PathLike[str], str, dict[str, DataFile], None],
        # filter_kwargs: Optional[dict[str, Any]] = None,
        # correction_kwargs: Optional[dict[str, Any]] = None,
    ) -> None:

        self.hdx_spec = hdx_spec
        self.data_files: dict[str, DataFile] = {}

        if isinstance(data_src, (os.PathLike, str)):
            data_src = Path(data_src) or Path(".")
            for name, spec in self.hdx_spec['data_files'].items():
                datafile = DataFile(name=name,
                                    filepath_or_buffer= data_src / spec['filename'],
                                    **{k: v for k, v in spec.items() if k != 'filename'},
                                    )
                self.data_files[name] = datafile

        elif isinstance(data_src, dict):
            self.data_files = data_src
        else:
            raise TypeError(
                f"Invalid data type {type(data_src)!r}, must be path or dict"
            )

    def load_hdxmset(self) -> HDXMeasurementSet:
        hdxm_list = [self.load_hdxm(state) for state in self.hdx_spec["states"].keys()]
        return HDXMeasurementSet(hdxm_list)

    def load_peptides(self, state: Union[str, int], peptides: str) -> pd.DataFrame:
        state = self.states[state] if isinstance(state, int) else state
        peptide_spec = self.hdx_spec["states"][state]['peptides'][peptides]
        filter_fields = {'state', 'exposure', 'query', 'dropna'}

        df = self.data_files[peptide_spec['data_file']].data
        peptides = filter_peptides(df, **{k: v for k, v in peptide_spec.items() if k in filter_fields})

        return peptides

    def load_hdxm(self, state: Union[str, int]) -> HDXMeasurement:
        state = self.states[state] if isinstance(state, int) else state
        peptide_spec = self.hdx_spec["states"][state]['peptides']
        metadata = self.hdx_spec["states"][state]['metadata']

        peptides = self.load_peptides(state, 'experiment')
        fd_peptides = self.load_peptides(state, 'FD_control') if 'FD_control' in peptide_spec else None
        nd_peptides = self.load_peptides(state, 'ND_control') if 'ND_control' in peptide_spec else None

        if fd_peptides is not None and 'be_percent' in metadata:
            peptides = correct_d_uptake(peptides)
            back_exchange = metadata['be_percent'] / 100.
            peptides["rfu"] = peptides["uptake"] / ((1 - back_exchange) * peptides['ex_residues'])
            peptides["uptake_corrected"] = peptides["uptake"] / (1 - back_exchange)
        else:
            peptides = apply_control(peptides, fd_peptides, nd_peptides)
            peptides = correct_d_uptake(peptides, drop_first=cfg.analysis.drop_first, )

        hdxm = HDXMeasurement(peptides, name=state, **metadata)

        return hdxm

    @property
    def correction_kwargs(self):
        kwargs = {
            "drop_first": cfg.analysis.drop_first,
            "d_percentage": self.hdx_spec['metadata'].get("d_percentage", 100.)
        }

        # todo:
        # if 'corrections' in self.hdx_spec:
        # ...

        return kwargs

    @property
    def states(self) -> list[str]:
        return list(self.hdx_spec['states'].keys())