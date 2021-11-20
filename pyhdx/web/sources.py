import pandas as pd
import param

from pyhdx import TorchFitResult
from pyhdx.fitting import RatesFitResult
from pyhdx.models import HDXMeasurement, HDXMeasurementSet


#todo refactor module to models?


class AppSourceBase(param.Parameterized):
    """Base class for sources"""

    _type = 'base'

    updated = param.Event()

    def get(self):
        raise NotImplementedError()


class TableSource(AppSourceBase):

    tables = param.Dict({})

    _type = 'table'

    def get(self):
        if len(self.tables) == 1:
            return next(iter(self.tables.values))
        else:
            raise ValueError("TableSource has multiple tables, use `get_table`")

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

    _type = 'pyhdx'

    # table options are (table_name, (opts)):  (General: <quantity>_<specifier> -> opts[qty] for colors
    # peptides
    # rfu_residues (rfu)
    # rates
    # dG_fits (dG)
    # ddG_comparison

    #data objects: (perhaps these should be on the main controller? to make the source more general / serialiable?
    hdxm_objects = param.Dict({})
    rate_results = param.Dict({})  # dataframes?
    dG_fits = param.Dict({})
    #cmaps = param.Dict(default=CMAP_DEFAULTS)

    def from_file(self):
        pass
        # todo load hdxms first
        #then use those to reload dG results

    def add(self, obj, name):  # todo Name is None and use obj name?
        if isinstance(obj, HDXMeasurement):
            self.hdxm_objects[name] = obj
            self.param.trigger('hdxm_objects')
        elif isinstance(obj, TorchFitResult):
            self.dG_fits[name] = obj
            self.param.trigger('dG_fits')
        elif isinstance(obj, RatesFitResult):
            self.rate_results[name] = obj
            self.param.trigger('rate_results')
        else:
            raise ValueError(f"Unsupported object {obj!r}")

    @property
    def hdx_set(self):
        return HDXMeasurementSet(list(self.hdxm_objects.values()))

    @param.depends('hdxm_objects', watch=True)
    def _hdxm_objects_updated(self):
        combined = pd.concat([hdxm.data for hdxm in self.hdxm_objects.values()], axis=1,
                             keys=self.hdxm_objects.keys(), names=['state', 'quantity'])  #todo 'state' or 'name' or 'protein_state'?
        # todo catch valueerror duplicate entries
        pivoted = combined \
            .stack(level=0) \
            .pivot(index='id', columns=['state', 'exposure']) \
            .reorder_levels(['state', 'exposure', 'quantity'], axis=1) \
            .sort_index(axis=1)

        self.tables['peptides'] = pivoted  # level 3 multiindex

        #RFU per residue per exposure
        dfs = [hdxm.rfu_residues for hdxm in self.hdxm_objects.values()]
        combined = pd.concat(dfs, axis=1, keys=self.hdxm_objects.keys(), names=['state', 'exposure'])
        self.tables['rfu_residues'] = combined

        # todo this erorrs: self.param.trigger('tables')
        self.updated = True

    @param.depends('dG_fits', watch=True)
    def _fit_results_updated(self):  #todo method name / result dicts names
        combined = pd.concat([fit_result.output for fit_result in self.dG_fits.values()], axis=1,
                             keys=self.dG_fits.keys(), names=['fit_ID', 'state', 'quantity'])
        self.tables['dG_fits'] = combined

        self.updated = True

        # todo add d_exp etc
        #cached?:

    @param.depends('rate_results', watch=True)
    def _rates_results_updated(self):
        combined = pd.concat([fit_result.output for fit_result in self.rate_results.values()], axis=1,
                             keys=self.rate_results.keys(), names=['guess_ID', 'state', 'quantity'])
        self.tables['rates'] = combined

        self.updated = True


