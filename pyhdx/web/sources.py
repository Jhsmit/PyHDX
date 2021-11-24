import pandas as pd
import param
import numpy as np

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
            #self.dG_fits[name] = obj
            self._add_dG_fit(obj, name)
            #self.param.trigger('dG_fits')
        elif isinstance(obj, RatesFitResult):
            self.rate_results[name] = obj
            self.param.trigger('rate_results')
        else:
            raise ValueError(f"Unsupported object {obj!r}")

    @property
    def hdx_set(self):
        return HDXMeasurementSet(list(self.hdxm_objects.values()))

    @param.depends('hdxm_objects', watch=True)
    def _hdxm_objects_updated(self):  # todo change to hdxm object added?
        combined = pd.concat([hdxm.data for hdxm in self.hdxm_objects.values()], axis=1,
                             keys=self.hdxm_objects.keys(), names=['state', 'quantity'])  #todo 'state' or 'name' or 'protein_state'?
        # todo catch valueerror duplicate entries
        # todo this pivot reuses 'state' column entries which gives wrong 'state' in name in final index
        # should be user-entered state name
        # also: make index dtype int
        pivoted = combined \
            .stack(level=0) \
            .pivot(index='peptide_id', columns=['state', 'exposure']) \
            .reorder_levels(['state', 'exposure', 'quantity'], axis=1) \
            .sort_index(axis=1)

        self.tables['peptides'] = pivoted  # level 3 multiindex

        # RFU per residue per exposure
        dfs = [hdxm.rfu_residues for hdxm in self.hdxm_objects.values()]
        combined = pd.concat(dfs, axis=1, keys=self.hdxm_objects.keys(), names=['state', 'exposure'])
        self.tables['rfu_residues'] = combined

        # todo this erorrs: self.param.trigger('tables')
        self.updated = True

    def _add_dG_fit(self, fit_result, name):
        # Add deltaG values table (+ covariances etc)
        df = fit_result.output.copy()
        tuples = [(name, *tup) for tup in df.columns]
        columns = pd.MultiIndex.from_tuples(tuples, names=['fit_ID', 'state', 'quantity'])
        df.columns = columns

        if 'dG_fits' in self.tables:
            current = self.tables['dG_fits']
            new = pd.concat([current, df], axis=1)
        else:
            new = df
        self.tables['dG_fits'] = new

        # Add calculated d-uptake values (#todo add method on FitResults object that does this?)
        timepoints = fit_result.hdxm_set.timepoints
        tmin = np.log10(timepoints[np.nonzero(timepoints)].min())
        tmax = np.log10(timepoints.max())
        pad = 0.05 * (tmax - tmin)  # 5% padding percentage

        tvec = np.logspace(tmin - pad, tmax + pad, num=100, endpoint=True)
        d_calc = fit_result(tvec)

        # Reshape the c_calc numpy array (Ns x Np x Nt to pandas dataframe (index: Ns, columns: multiiindex Ns, Np)
        Ns, Np, Nt = d_calc.shape
        reshaped = d_calc.reshape(Ns * Np, Nt)
        columns = pd.MultiIndex.from_product(
            [[name], fit_result.hdxm_set.names, np.arange(Np), ['d_calc']],
            names=['Fit one', 'state', 'peptide_id', 'quantity'])
        index = pd.Index(tvec, name='exposure')
        df = pd.DataFrame(reshaped.T, index=index, columns=columns)
        df = df.loc[:, (df != 0).any(axis=0)]  # remove zero columns, replace with NaN when possible

        if 'd_calc' in self.tables:
            current = self.tables['d_calc']
            new = pd.concat([current, df], axis=1)
        else:
            new = df
        self.tables['d_calc'] = new

        # Add losses df
        df = fit_result.losses.copy()

        tuples = [(name, column) for column in df.columns]  # losses df is not multiindex
        columns = pd.MultiIndex.from_tuples(tuples, names=['fit_ID', 'loss_type'])
        df.columns = columns

        if 'loss' in self.tables:
            current = self.tables['loss']
            new = pd.concat([current, df], axis=1)
        else:
            new = df
        self.tables['loss'] = new

        #Add MSE per peptide df

        squared_errors = fit_result.get_squared_errors()

        dfs = {}
        for mse_sample, hdxm in zip(squared_errors, fit_result.hdxm_set):
            peptide_data = hdxm[0].data
            mse = np.mean(mse_sample, axis=1)
            # Indexing of mse_sum with Np to account for zero-padding
            data_dict = {'start': peptide_data['start'], 'end': peptide_data['end'], 'peptide_mse': mse[:hdxm.Np]}
            dfs[hdxm.name] = pd.DataFrame(data_dict)

        mse_df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1)
        mse_df.index.name = 'peptide_id'
        tuples = [(name, *tup) for tup in mse_df.columns]
        columns = pd.MultiIndex.from_tuples(tuples, names=['fit_ID', 'state', 'quantity'])
        mse_df.columns = columns

        self.tables['peptide_mse'] = mse_df

        self.dG_fits[name] = fit_result

        self.updated = True

#    @param.depends('dG_fits', watch=True)
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


