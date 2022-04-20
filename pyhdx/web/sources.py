import os
import urllib.request

import numpy as np
import pandas as pd
import param

from pyhdx import TorchFitResult, TorchFitResultSet
from pyhdx.fitting import RatesFitResult
from pyhdx.models import HDXMeasurement, HDXMeasurementSet
from pyhdx.support import multiindex_astype, multiindex_set_categories, hash_dataframe
from pyhdx.config import cfg


class Source(param.Parameterized):
    """Base class for sources"""

    _type = "base"

    updated = param.Event()

    def get(self):
        raise NotImplementedError()


class TableSource(Source):

    tables = param.Dict(default={}, doc="Dictionary of tables (pd.DataFrames)")

    hashes = param.Dict(default={}, doc="Dictionary of table hashes")

    _type = "table"

    def get(self):
        if len(self.tables) == 1:
            return next(iter(self.tables.values))
        else:
            raise ValueError("TableSource has multiple tables, use `get_table`")

    def add_table(self, table: str, df: pd.DataFrame):
        table_hash = hash_dataframe(df)
        self.hashes[table] = table_hash
        self.tables[table] = df

        # todo self.updated = True?

    def get_table(self, table):
        df = self.tables.get(table, None)

        return df

    def get_tables(self):
        """
        Returns the list of tables available on this source.
        Returns
        -------
        list
            The list of available tables on this source.
        """

        return list(self.tables.keys())


class PyHDXSource(TableSource):
    _type = "pyhdx"

    # see readme/tables_list for tables and their indexes

    hdxm_objects = param.Dict({})
    rate_results = param.Dict({})  # dict of rate fitting / guesses results
    dG_fits = param.Dict({})  # dict of torch fit result objects

    def from_file(self):
        pass
        # todo load hdxms first
        # then use those to reload dG results

    def add(self, obj, name):  # todo Name is None and use obj name?
        if isinstance(obj, HDXMeasurement):
            self._add_hdxm_object(obj, name)
        elif isinstance(obj, (TorchFitResult, TorchFitResultSet)):
            self._add_dG_fit(obj, name)
        elif isinstance(obj, RatesFitResult):
            self._add_rates_fit(obj, name)
        else:
            raise ValueError(f"Unsupported object {obj!r}")

    @property
    def hdx_set(self):
        return HDXMeasurementSet(list(self.hdxm_objects.values()))

    @property
    def names(self):
        """returns the names of all HDX Measurment objects loaded"""
        return list(self.hdxm_objects.keys())

    def _add_rates_fit(self, rates_result, name):
        df = rates_result.output.copy()

        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples, names=["guess_ID", "state", "quantity"]
        )
        df.columns = columns
        self._add_table(df, "rates")
        self.rate_results[name] = rates_result
        self.param.trigger("rate_results")
        self.updated = True

    def _add_hdxm_object(
        self, hdxm, name
    ):  # where name is new 'protein state' entry (or used for state (#todo clarify))
        # Add peptide data
        df = hdxm.data_wide.copy()
        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples, names=["state", "exposure", "quantity"]
        )
        df.columns = columns
        self._add_table(df, "peptides")

        # Add rfu per residue data
        df = hdxm.rfu_residues
        tuples = [(name, column, "rfu") for column in df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples, names=["state", "exposure", "quantity"]
        )
        df.columns = columns
        self._add_table(df, "rfu_residues")

        self.hdxm_objects[name] = hdxm
        self.param.trigger("hdxm_objects")  # protein controller listens here
        self.updated = True

    def _add_dG_fit(self, fit_result, name):
        # Add dG values table (+ covariances etc)
        df = fit_result.output.copy()
        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples, names=["fit_ID", "state", "quantity"]
        )
        df.columns = columns
        self._add_table(df, "dG_fits")

        # Add calculated d-uptake values
        df = fit_result.get_dcalc()

        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples, names=["fit_ID", "state", "peptide_id", "quantity"],
        )
        df.columns = columns

        self._add_table(df, "d_calc")

        # Add losses df
        df = fit_result.losses.copy()
        if df.columns.nlevels == 1:
            tuples = [(name, column) for column in df.columns]
            columns = pd.MultiIndex.from_tuples(tuples, names=["fit_ID", "loss_type"])
        else:
            tuples = [(name, *tup) for tup in df.columns]
            columns = pd.MultiIndex.from_tuples(
                tuples, names=["fit_ID", "state", "loss_type"]
            )

        df.columns = columns
        self._add_table(df, "loss")

        # Add MSE per peptide df
        # current bug: convert dtypes drop column names: https://github.com/pandas-dev/pandas/issues/41435
        # use before assigning column names
        mse_df = fit_result.get_peptide_mse().convert_dtypes()
        # mse_df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1).convert_dtypes()
        mse_df.index.name = "peptide_id"
        tuples = [(name, *tup) for tup in mse_df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples, names=["fit_ID", "state", "quantity"]
        )
        mse_df.columns = columns
        self._add_table(mse_df, "peptide_mse")

        self.dG_fits[name] = fit_result
        self.updated = True

    def _add_table(
        self, df, table, categorical=True
    ):  # TODO add_table is (name, dataframe)
        """

        :param df:
        :param table: name of the table
        :param categorical: True if top level of multiindex should be categorical
        :return:
        """

        if table in self.tables:
            current = self.tables[table]
            new = pd.concat([current, df], axis=1, sort=True)
            categories = list(current.columns.unique(level=0)) + list(
                df.columns.unique(level=0)
            )
        else:
            new = df
            categories = list(df.columns.unique(level=0))
        if categorical:
            new.columns = multiindex_astype(new.columns, 0, "category")
            new.columns = multiindex_set_categories(
                new.columns, 0, categories, ordered=True
            )

        self.add_table(table, new)


class PDBSource(Source):

    _type = "pdb"

    pdb_files = param.Dict({}, doc="Dictionary with id: pdb_string pdb file entries")

    hashes = param.Dict({})

    max_entries = param.Number(
        1,
        doc="set maximum size for pdb files. set to none for infinite size. set to one for single pdb mode",
    )

    def add_from_pdb(self, pdb_id):
        self._make_room()
        url = f"http://files.rcsb.org/download/{pdb_id}.pdb"
        # with urllib.request.urlopen(url) as response:
        #     pdb_string = response.read().decode()

        self.pdb_files[pdb_id] = url
        self.hashes[pdb_id] = hash(url)
        self.updated = True

    def add_from_string(self, pdb_string, pdb_id):
        """Adds a PDB file from a string"""
        self._make_room()
        url = f"assets/{pdb_id}.pdb"

        (cfg.assets_dir / f"{pdb_id}.pdb").write_text(pdb_string)

        self.pdb_files[pdb_id] = url
        self.hashes[pdb_id] = hash(url)
        self.updated = True

    def _make_room(self):
        """removes first entry of pdb_files dict if its at max capacity"""
        if self.max_entries is None:
            pass
        elif len(self.pdb_files) == self.max_entries:
            key = next(iter(self.pdb_files))
            del self.pdb_files[key]
            del self.hashes[key]
            pdb_pth = cfg.assets_dir / f"{key}.pdb"
            # if pdb_pth.exists():
            #     os.remove(pdb_pth)

    def get(self):
        """returns the first entry in the pdb source"""
        return next(iter(self.pdb_files.values()))

    def get_pdb(self, pdb_id):
        return self.pdb_files[pdb_id]
