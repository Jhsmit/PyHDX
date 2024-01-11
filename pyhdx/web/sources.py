import json
import uuid
from collections import defaultdict
from typing import Optional

import pandas as pd
import param

from pyhdx import TorchFitResult, TorchFitResultSet
from pyhdx.fitting import RatesFitResult, DUptakeFitResultSet
from pyhdx.models import HDXMeasurement, HDXMeasurementSet
from pyhdx.support import multiindex_astype, multiindex_set_categories, hash_dataframe
from pyhdx.config import cfg

# <table_name>: {'cmap_field': <table_column_name>, cmap_opt: <cmap_opt_name>
from pyhdx.web.utils import fix_multiindex_dtypes

TABLE_INFO = {
    "rfu": {"cmap_field": "rfu", "cmap_opt": "rfu_cmap"},
    "d_uptake": {"cmap_field": "d_uptake", "cmap_opt": "d_uptake_cmap"},
    "dd_uptake": {"cmap_field": "dd_uptake", "cmap_opt": "dd_uptake_cmap"},
    "dG": {"cmap_field": "dG", "cmap_opt": "dG_cmap"},
    "drfu": {"cmap_field": "drfu", "cmap_opt": "drfu_cmap"},
    "ddG": {"cmap_field": "ddG", "cmap_opt": "ddG_cmap"},
    #'rfu': {'cmap_field': 'rfu', 'cmap_opt': 'rfu_cmap'},
}


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
        if len(self.tables) == 0:
            return None
        elif len(self.tables) == 1:
            return next(iter(self.tables.values()))

        else:
            raise ValueError("TableSource has multiple tables, use `get_table`")

    def set(self, df: pd.DataFrame) -> None:
        if len(self.tables) > 1:
            raise ValueError(
                "Can only use the `set` method when the number of tables is below 1. "
                "Use `add_table`"
            )

        table = next(iter(self.tables.keys())) if self.tables else "main"
        self.add_table(table, df)

    def add_table(self, table: str, df: pd.DataFrame) -> None:
        table_hash = hash_dataframe(df)
        self.hashes[table] = table_hash
        self.tables[table] = df

        # todo self.updated = True? (yes because right now this is done manually (?))
        # although would be good to have an option to not trigger (or context manager)
        # when adding multiple tables in batch and not wanted to update

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
    d_uptake_results = param.Dict({})
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
        elif isinstance(obj, DUptakeFitResultSet):
            self._add_duptake_fit(obj, name)
        else:
            raise ValueError(f"Unsupported object {obj!r}")

    @property
    def hdx_set(self):
        return HDXMeasurementSet(list(self.hdxm_objects.values()))

    @property
    def names(self):
        """returns the names of all HDX Measurment objects loaded"""
        return list(self.hdxm_objects.keys())

    def _add_duptake_fit(self, d_uptake_result, name):
        df = d_uptake_result.output
        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples, names=["D_uptake_fit_ID", "state", "exposure", "quantity"]
        )

        df.columns = fix_multiindex_dtypes(columns)
        self._add_table(df, "d_uptake")
        self.d_uptake_results[name] = d_uptake_result
        self.param.trigger("d_uptake_results")  # todo no listeners probably
        self.updated = True

    def _add_rates_fit(self, rates_result, name):
        df = rates_result.output.copy()

        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(tuples, names=["guess_ID", "state", "quantity"])
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
        columns = pd.MultiIndex.from_tuples(tuples, names=["state", "exposure", "quantity"])
        df.columns = columns
        self._add_table(df, "peptides")

        # Add rfu per residue data
        # todo perhaps this combined df should be directly supplied by `hdxm`
        rfu = hdxm.rfu_residues
        columns = pd.MultiIndex.from_tuples(
            [(name, col, "rfu") for col in rfu.columns],
            names=["state", "exposure", "quantity"],
        )
        rfu.columns = columns

        rfu_sd = hdxm.rfu_residues_sd
        columns = pd.MultiIndex.from_tuples(
            [(name, col, "rfu_sd") for col in rfu_sd.columns],
            names=["state", "exposure", "quantity"],
        )
        rfu_sd.columns = columns

        combined = pd.concat([rfu, rfu_sd], axis=1).sort_index(axis=1)

        self._add_table(combined, "rfu")

        self.hdxm_objects[name] = hdxm
        self.param.trigger("hdxm_objects")  # protein controller listens here
        self.updated = True

    def _add_dG_fit(self, fit_result, name):
        # Add dG values table (+ covariances etc)
        df = fit_result.output.copy()
        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(tuples, names=["fit_ID", "state", "quantity"])
        df.columns = columns
        self._add_table(df, "dG")

        # Add calculated d-uptake values
        df = fit_result.get_dcalc()

        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(
            tuples,
            names=["fit_ID", "state", "peptide_id", "quantity"],
        )
        df.columns = columns

        self._add_table(df, "d_calc")

        # Add losses df
        df = fit_result.losses.copy()
        if df.columns.nlevels == 1:
            tuples = [(name, "*", column) for column in df.columns]
            columns = pd.MultiIndex.from_tuples(tuples, names=["fit_ID", "state", "loss_type"])
        else:
            tuples = [(name, *tup) for tup in df.columns]
            columns = pd.MultiIndex.from_tuples(tuples, names=["fit_ID", "state", "loss_type"])

        df.columns = columns
        self._add_table(df, "loss")

        # Add MSE per peptide df
        # current bug: convert dtypes drop column names: https://github.com/pandas-dev/pandas/issues/41435
        # use before assigning column names
        mse_df = fit_result.get_peptide_mse().convert_dtypes()
        # mse_df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1).convert_dtypes()
        mse_df.index.name = "peptide_id"
        tuples = [(name, *tup) for tup in mse_df.columns]
        columns = pd.MultiIndex.from_tuples(tuples, names=["fit_ID", "state", "quantity"])
        mse_df.columns = columns
        self._add_table(mse_df, "peptide_mse")

        self.dG_fits[name] = fit_result
        self.updated = True

    def _add_table(self, df, table, categorical=True):  # TODO add_table is (name, dataframe)
        """

        :param df:
        :param table: name of the table
        :param categorical: True if top level of multiindex should be categorical
        :return:
        """

        if table in self.tables:
            current = self.tables[table]
            new = pd.concat([current, df], axis=1, sort=True)
            categories = list(current.columns.unique(level=0)) + list(df.columns.unique(level=0))
        else:
            new = df
            categories = list(df.columns.unique(level=0))
        if categorical:
            new.columns = multiindex_astype(new.columns, 0, "category")
            new.columns = multiindex_set_categories(new.columns, 0, categories, ordered=True)

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


class DictSource(Source):
    """Source for (metadata) dictionaries"""

    _type = "dict"

    _items = param.Dict({})

    hashes = param.Dict({})

    def set(self, item: dict, name: Optional[str] = None):
        if not isinstance(item, dict):
            raise TypeError(f"Invalid type of 'item', must be {dict!r}, got {type(item)!r}")
        # self.make_room()
        name = name or f"_item_{uuid.uuid4()}"  # self.new_key()
        self._items[name] = item
        self.hashes[name] = self.hash_item(item)

    # def set_value(self, name: str, key: Any, value: Any,):
    #     d = self._items[name]
    #
    #
    #     d[key] = value
    #     self.hashes[name] = self.hash_item(self._items[name])

    def hash_item(self, item: dict) -> int:
        return hash(json.dumps(item))

    def update(self) -> None:
        self.hashes = {key: self.hash_item(item) for key, item in self._items}
        self.updated = True

    # todo does not match base object
    def get(self, name: str) -> Optional[dict]:
        if name not in self._items:
            item = defaultdict(dict)
            self._items[name] = item
            return item

        name = name or next(iter(self.keys()))
        return self._items[name]
