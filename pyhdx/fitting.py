from pyhdx.support import get_reduced_blocks, temporary_seed
from pyhdx.models import Protein
from pyhdx.fitting_torch import DeltaGFit, TorchSingleFitResult, TorchBatchFitResult
from pyhdx.fit_models import SingleKineticModel, OneComponentAssociationModel, TwoComponentAssociationModel, OneComponentDissociationModel, \
    TwoComponentDissociationModel
from scipy import constants
from scipy.optimize import fsolve
import torch
import numpy as np
from symfit import Fit, Variable, Parameter, exp, Model, CallableModel
from symfit.core.minimizers import DifferentialEvolution, Powell
from collections import namedtuple
from functools import reduce, partial
from operator import add
from dask.distributed import Client
import warnings


EmptyResult = namedtuple('EmptyResult', ['chi_squared', 'params'])
er = EmptyResult(np.nan, {k: np.nan for k in ['tau1', 'tau2', 'r']})


def fit_kinetics(t, d, model, chisq_thd):
    """
    Fit time kinetics with two time components and corresponding relative amplitude.

    Parameters
    ----------
    t : :class:`~numpy.ndarray`
        Array of time points
    d : :class:`~numpy.ndarray`
        Array of uptake values
    chisq_thd : :obj:`float`
        Threshold chi squared above which the fitting is repeated with the Differential Evolution algorithm.

    Returns
    -------
    res : :class:`~symfit.FitResults`
        Symfit fitresults object.
    """
    if np.any(np.isnan(d)):
        raise ValueError('There shouldnt be NaNs anymore')
        er = EmptyResult(np.nan, {p.name: np.nan for p in model.sf_model.params})
        return er

    model.initial_guess(t, d)
    with temporary_seed(43):
        fit = Fit(model.sf_model, t, d, minimizer=Powell)
        res = fit.execute()

        if not check_bounds(res) or np.any(np.isnan(list(res.params.values()))) or res.chi_squared > chisq_thd:
            fit = Fit(model.sf_model, t, d, minimizer=DifferentialEvolution)
            #grid = model.initial_grid(t, d, step=5)
            res = fit.execute()

    return res


def check_bounds(fit_result):
    for param in fit_result.model.params:
        value = fit_result.params[param.name]
        if value < param.min:
            return False
        elif value > param.max:
            return False
    return True


def fit_global(data, model):
    fit = Fit(model.sf_model, **data)
    res = fit.execute()
    return res


# Defaults for PyTorch optimizations
optimizer_defaults = {
    'SGD': {
        'lr': 10,
        'momentum': 0.5,
        'nesterov': True
    },
}


def run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs, model, criterion, regularizer,
                  epochs=100000, patience=50, stop_loss=0.05):

    optimizer_obj = optimizer_klass(model.parameters(), **optimizer_kwargs)

    np.random.seed(43)
    torch.manual_seed(43)

    mse_loss_list = [np.inf]
    total_loss_list = [np.inf]

    def closure():
        output = model(*inputs)
        loss = criterion(output, output_data)
        mse_loss_list.append(loss.detach())
        reg_loss = regularizer(model.deltaG)
        total_loss = loss + reg_loss
        total_loss_list.append(total_loss.detach())
        total_loss.backward()
        return total_loss

    stop = 0
    for epoch in range(epochs):
        optimizer_obj.zero_grad()
        loss = optimizer_obj.step(closure)

        diff = total_loss_list[-2] - total_loss_list[-1]
        if diff < stop_loss:
            stop += 1
            if stop > patience:
                break
        else:
            stop = 0

    #par = model.deltaG.detach().numpy()
    return np.array(mse_loss_list), np.array(total_loss_list), model


class KineticsFitting(object):

    def __init__(self, series, bounds=None, temperature=None, pH=None, c_term=None, cluster=None):
        #todo perhaps this whole object should dissapear in favour of series
        # and make fitting functional instead of object oriented
        self.series = series
        self.bounds = bounds or self._get_bounds()

        #todo temperature, pH, c_term from series?
        self.temperature = temperature
        self.pH = pH
        self.c_term = c_term
        self.cluster = cluster

    def set_k_int(self):
        """
        Sets intrinsic rates of exchange on the KineticsSeries protein object

        Returns
        -------

        Nothing

        """

        self.series.coverage.protein.set_k_int(self.temperature, self.pH)

    def _get_bounds(self):
        #todo document
        #todo this is now causing confusion when getting protection factors as output
        times = self.series.timepoints
        nonzero_times = times[np.nonzero(times)]
        t_first = np.min(nonzero_times)
        t_last = np.max(nonzero_times)
        b_upper = 50*np.log(2) / t_first
        b_lower = np.log(2) / (t_last * 50)

        return b_lower, b_upper

    def _prepare_wt_avg_fit(self, model_type='association'):
        models = []
        intervals = []  # Intervals; (start, end); (inclusive, exclusive)
        d_list = []
        series = self.series

        arr = series.scores_stack.T
        i = 0
        # because intervals are inclusive, exclusive we need to add an extra entry to r_number for the final exclusive bound
        r_excl = np.append(series.coverage.r_number, [series.coverage.r_number[-1] + 1])

        for bl in series.coverage.block_length:
            d = arr[i]
            if np.all(np.isnan(d)):  # Skip non-coverage blocks
                i += bl
                continue
            intervals.append((r_excl[i], r_excl[i + bl]))
            d_list.append(d)
            if model_type == 'association':
                model = TwoComponentAssociationModel(self.bounds)
            elif model_type == 'dissociation':
                model = TwoComponentDissociationModel(self.bounds)
            else:
                raise ValueError('Invalid model type {}'.format(model_type))
            # res = EmptyResult(np.nan, {p.name: np.nan for p in model.sf_model.params})

            models.append(model)
            i += bl  # increment in block length does not equal move to the next start position

        return d_list, intervals, models

    async def weighted_avg_fit_async(self, chisq_thd=20, model_type='association', pbar=None):
        """
        Block length _should_ be equal to the block length of all measurements in the series, provided that all coverage
        is the same

        Parameters
        ----------
        chisq_max

        Returns
        -------

        """

        d_list, intervals, models = self._prepare_wt_avg_fit(model_type=model_type)
        fit_func = partial(fit_kinetics, self.series.timepoints)
        client = await Client(self.cluster)
        futures = client.map(fit_func, d_list, models, chisq_thd=chisq_thd)
        if pbar:
            pbar.num_tasks = len(d_list)  #this call assignment might also need unlocked()
            await pbar.run(futures)

        results = client.gather(futures)

        fit_result = KineticsFitResult(self.series, intervals, results, models)
        return fit_result

    def weighted_avg_fit(self, chisq_thd=20, model_type='association', pbar=None, callbacks=None):
        """
        Block length _should_ be equal to the block length of all measurements in the series, provided that all coverage
        is the same

        Parameters
        ----------
        chisq_max

        Returns
        -------

        """

        d_list, intervals, models = self._prepare_wt_avg_fit(model_type=model_type)
        if pbar:
            self.num_tasks = len(d_list)
            inc = pbar.increment
        else:
            inc = lambda: None

        results = []
        for d, model in zip(d_list, models):
            result = fit_kinetics(self.series.timepoints, d, model, chisq_thd=chisq_thd)
            inc()
            results.append(result)

        fit_result = KineticsFitResult(self.series, intervals, results, models)
        return fit_result

    def guess_deltaG(self, guess_rates):
        #todo do some checks on the index of supplied guess_rates
        protein = self.series.coverage.protein
        p_guess = (protein['k_int'] / guess_rates['rate']) - 1
        p_guess.clip(0., None, inplace=True)  # Some initial guesses will have negative PF values
        deltaG = np.log(p_guess) * constants.R * self.temperature

        bools = ~np.isfinite(deltaG)
        idx = np.where(np.diff(bools))[0]
        i = 0 if np.isfinite(deltaG.iloc[0]) else 1  # Determine if guesses start with coverage/data or not
        for start, stop in zip(idx[i::2], idx[1 + i::2]):
            replacement = np.linspace(deltaG.iloc[start], deltaG.iloc[stop + 1], endpoint=True,
                                      num=stop - start + 2)
            deltaG.iloc[start + 1: stop + 1] = replacement[1:-1]

        # Guesses end with NaN block:
        if (i + len(idx)) % 2 == 1:
            deltaG.iloc[idx[-1]:] = deltaG.iloc[idx[-1]]

        return deltaG

    def _initial_guess(self, initial_guess, protein):
        #todo refactor prepare_initial_guess?
        #todo deprecate
        """

        Parameters
        ----------
        initial_rates: :class:`pyhdx.models.Protein`
            Must have rate column
        protein :class:`pyhdx.models.Protein`
            has r_number, k_int as keys, values are numpy arrays

        Returns
        -------

        """

        p_guess = (protein['k_int'] / initial_guess['rate']) - 1

        # Replace NaN (no coverage / prolines) with logspace interpolation between edges
        bools = np.logical_or(np.isnan(p_guess), p_guess == 0.)
        idx = np.where(np.diff(bools))[0]
        i = 1 if np.isnan(p_guess.iloc[0]) else 0
        for start, stop in zip(idx[i::2], idx[1 + i::2]):
            replacement = np.logspace(np.log10(p_guess.iloc[start]), np.log10(p_guess.iloc[stop + 1]), endpoint=True,
                                      num=stop - start + 2)
            p_guess.iloc[start + 1:stop + 1] = replacement[1:-1]

        return p_guess

    #todo might make more sense to have initial result as deltaG vecotor as input
    def global_fit(self, initial_result, r1=2, epochs=100000, patience=50, stop_loss=0.05,
                   optimizer='SGD', **optimizer_kwargs):
        #todo @tejas: Missing docstring
        """Pytorch global fitting"""

        deltaG_par, inputs, output_data = self.setup_global_fit(initial_result)
        model = DeltaGFit(deltaG_par)
        criterion = torch.nn.MSELoss(reduction='sum')

        # Take default optimizer kwargs and update them with supplied kwargs
        optimizer_kwargs = {**optimizer_defaults.get(optimizer, {}), **optimizer_kwargs}  # Take defaults and override with user-specified
        optimizer_klass = getattr(torch.optim, optimizer)

        def regularizer(param):
            return r1 * torch.mean(torch.abs(param[:-1] - param[1:]))

        # returned_model is the same object as model
        mse_loss, total_loss, returned_model = run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs,
                                                             model, criterion, regularizer, epochs=epochs,
                                                             patience=patience, stop_loss=stop_loss)

        result = TorchSingleFitResult(self, model,
                                      mse_loss=mse_loss, total_loss=total_loss)

        return result

    def setup_global_fit(self, initial_result, r1=2):
        #todo @tejas: Missing docstring
        """Pytorch global fitting"""

        if 'k_int' not in self.series.coverage.protein:
            self.set_k_int()

        dtype = torch.float64

        # Prepare input data in the correct shapes for fitting
        temperature = torch.tensor([self.temperature], dtype=dtype)
        X = torch.tensor(self.series.coverage.X, dtype=dtype) # Np x Nr
        k_int = torch.tensor(self.series.coverage['k_int'].to_numpy(), dtype=dtype).unsqueeze(-1)  # Nr x 1
        timepoints = torch.tensor(self.series.timepoints, dtype=dtype).unsqueeze(0)  # 1 x Nt
        inputs = [temperature, X, k_int, timepoints]

        # Prepare output data in the correct shape for fitting
        output_data = torch.tensor(self.series.uptake_corrected.T, dtype=dtype)

        # Get initial guess values for deltaG
        gibbs_values = self.series.coverage.apply_interval(self.guess_deltaG(initial_result)).to_numpy()
        if np.any(np.isnan(gibbs_values)):
            raise ValueError('NaN values in initial guess values')

        deltaG_par = torch.nn.Parameter(torch.Tensor(gibbs_values).unsqueeze(-1))

        return deltaG_par, inputs, output_data

    async def global_fit_async(self, initial_result, r1=2, epochs=100000, patience=50, stop_loss=0.05,
                               optimizer='SGD', **optimizer_kwargs):

        deltaG_par, inputs, output_data = self.setup_global_fit(initial_result)

        model = DeltaGFit(deltaG_par)
        criterion = torch.nn.MSELoss(reduction='sum')

        #todo repeated code with global_fit

        # Take default optimizer kwargs and update them with supplied kwargs
        optimizer_kwargs = {**optimizer_defaults.get(optimizer, {}), **optimizer_kwargs}  # Take defaults and override with user-specified
        optimizer_klass = getattr(torch.optim, optimizer)

        def regularizer(param):
            return r1 * torch.mean(torch.abs(param[:-1] - param[1:]))

        fit_func = partial(run_optimizer, inputs, output_data, optimizer_klass, optimizer_kwargs, model, criterion, regularizer,
                           epochs=epochs, patience=patience, stop_loss=stop_loss)
        client = await Client(self.cluster, asynchronous=True)

        future = client.submit(fit_func)

        # Get returned_model from Dask client which has updated params which are the fitted parameters
        mse_loss, total_loss, returned_model = await future

        await client.close()
        result = TorchSingleFitResult(self, returned_model, mse_loss=mse_loss, total_loss=total_loss)

        return result

    def weighted_avg_t50(self):
        """
        Calculates exchange rates based on weighted averaging followed by interpolation to determine half-time, which is
        then calculated to rates.

        Returns
        -------

        output: :class:`~np.ndarray`
            array with fields r_number, rate

        """
        #todo this is uing the soon to be depcrecated coverage object
        interpolated = np.array([np.interp(50, d_uptake, self.series.timepoints) for d_uptake in self.series.scores_stack.T])

        output = np.empty_like(interpolated, dtype=[('r_number', int), ('rate', float)])
        output['r_number'] = self.series.coverage.r_number
        output['rate'] = np.log(2) / interpolated

        protein = Protein(output, index='r_number')
        t50FitResult = namedtuple('t50FitResult', ['output'])

        result = t50FitResult(output=protein)

        return result

    def weighted_avg_linearize(self):
        rates = []
        output = np.empty_like(self.series.coverage.r_number, dtype=[('r_number', int), ('rate', float)])
        output['r_number'] = self.series.coverage.r_number

        for i, dpts in enumerate(self.series.scores_stack.T):
            if np.any(np.isnan(dpts)):
                output['rate'][i] = np.nan
                rates.append(np.nan)
            else:
                y_lin = np.log(1 - (dpts / 100))
                b = ~np.isnan(y_lin)
                try:
                    rate, offset = np.polyfit(self.series.timepoints[b], -y_lin[b], 1)
                except np.linalg.LinAlgError:
                    t50 = np.interp(50, dpts, self.series.timepoints)
                    rate = np.log(2) / t50
                output['rate'][i] = rate

        return output


class KineticsFitResult(object):
    """
    this fit results is only for wt avg fitting
    """
    def __init__(self, series, intervals, results, models):
        """
        each model with corresponding interval covers a region in the protein corresponding to r_number
        """
        assert len(results) == len(models)
#        assert len(models) == len(block_length)
        self.series = series
        self.r_number = series.coverage.r_number
        self.intervals = intervals  #inclusive, excluive
        self.results = results
        self.models = models

    @property
    def model_type(self):
        if np.all([isinstance(m, SingleKineticModel) for m in self.models]):
            return 'Single'
        else:
            raise ValueError('Unsupported model types')

    def __call__(self, timepoints):
        """call the result with timepoints to get fitted uptake per peptide back"""
        d_list = []
        if self.model_type == 'Single':
            for time in timepoints:
                p = self.get_p(time)
                p = np.nan_to_num(p)
                d = self.series.coverage.X.dot(p)
                d_list.append(d)
        elif self.model_type == 'Global':
            for time in timepoints:
                d = self.get_d(time)
                d_list.append(d)

        uptake = np.vstack(d_list).T
        return uptake

    def get_p(self, t):
        """
        Calculate P at timepoint t. Only for wt average type fitting results

        """
        p = np.full_like(self.r_number, fill_value=np.nan, dtype=float)
        for (s, e), result, model in zip(self.intervals, self.results, self.models):
            i0, i1 = np.searchsorted(self.r_number, [s, e])
            p[i0:i1] = model(t, **result.params)

        return p

    def get_d(self, t):
        "calculate d at timepoint t only for lsqkinetics (refactor glocal) type fitting results (scores per peptide)"
        d_arr = np.array([])
        for result, model in zip(self.results, self.models):
            d = np.array(model(t, **result.params)).flatten()
            d_arr = np.append(d_arr, d)
        return d_arr

    def __len__(self):
        return len(self.results)

    def __iter__(self):
        raise DeprecationWarning()
        iterable = [(r, m, b) for r, m, b in zip(self.results, self.models, self.block_length)]
        return iterable.__iter__()

    def get_param(self, name):
        """
        Get an array of parameter with name `name` from the fit result. The length of the array is equal to the
        number of amino acids.

        Parameters
        ----------
        name : :obj:`str`
            Name of the parameter to extract

        Returns
        -------
        par_arr : :class:`~numpy.ndarray`
            Array with parameter values

        """

        output = np.full_like(self.r_number, np.nan, dtype=float)
        for (s, e), result, model in zip(self.intervals, self.results, self.models):
            try:
                dummy_name = model.names[name]  ## dummy parameter name
                value = result.params[dummy_name]  #value is scalar
            except KeyError:   # is a lsqmodel funky town
                value = model.get_param_values(name, **result.params)  #value is vector
                #values = model.get_parameter(name)            #todo unify nomenclature param/parameter

            i0, i1 = np.searchsorted(self.r_number, [s, e])
            output[i0:i1] = value

        return output

    @property
    def rate(self):
        """Returns an array with the exchange rates"""
        output = np.full_like(self.r_number, np.nan, dtype=float)
        for (s, e), result, model in zip(self.intervals, self.results, self.models):
            rate = model.get_rate(**result.params)
            i0, i1 = np.searchsorted(self.r_number, [s, e])
            output[i0:i1] = rate
        return output

    @property
    def tau(self):
        """Returns an array with the exchange rates"""
        return 1 / self.rate

    def get_output(self, names):
        # change to property which gives all parameters as output
        dtype = [('r_number', int)] + [(name, float) for name in names]
        array = np.full_like(self.r_number, np.nan, dtype=dtype)
        array['r_number'] = self.r_number
        for name in names:
            try:
                array[name] = getattr(self, name)
            except AttributeError:
                array[name] = self.get_param(name)
        return array

    @property
    def output(self):
        array = self.get_output(['rate', 'k1', 'k2', 'r'])
        return Protein(array, index='r_number')


class BatchFitting(object):
    """Fit multiple datasets simultanuously in batch"""

    def __init__(self, states, guesses=None, cluster=None):
        #todo guesses as deltaG
        self.states = states

        #todo create Coverage object for the 3d case
        intervals = np.array([kf.series.coverage.interval for kf in self.states])
        self.interval = (intervals[:, 0].min(), intervals[:, 1].max())
        r_number = np.arange(*self.interval)
        self.r_number = r_number

        self.Ns = len(self.states)
        self.Nr = len(r_number)
        self.Np = np.max([kf.series.coverage.X.shape[0] for kf in self.states])
        self.Nt = np.max([len(kf.series.timepoints) for kf in self.states])

        self.guesses = guesses
        self.cluster = cluster

    def setup_fit(self):
        assert self.guesses is not None, 'Guesses are required to set up the fit'
        # Create numpy arrays with correct shapes as input data
        X = np.zeros((self.Ns, self.Np, self.Nr))
        D = np.zeros((self.Ns, self.Np, self.Nt))
        k_int = np.zeros((self.Ns, self.Nr))
        gibbs = np.full((self.Ns, self.Nr), fill_value=np.nan)  #todo default value for gibbs
        timepoints = np.zeros((self.Ns, self.Nt))

        # Set values for numpy input data
        for i, kf in enumerate(self.states):
            if 'k_int' not in kf.series.coverage.protein:
                kf.set_k_int()

            interval_sample = kf.series.coverage.interval
            # Indices of residues
            i0 = interval_sample[0] - self.interval[0]
            i1 = interval_sample[1] - self.interval[0]

            Npi = kf.series.coverage.X.shape[0]  # number of peptides in this particular state
            Nti = len(kf.series.timepoints) # number of timepoints in this particular state

            temperature = np.array([kf.temperature for kf in self.states])

            k_int_values = kf.series.coverage['k_int'].to_numpy()
            k_int[i, i0:i1] = k_int_values

            np.zeros((self.Ns, self.Nr))

            gibbs_values = kf.series.coverage.apply_interval(kf.guess_deltaG(self.guesses[i])).to_numpy()
            gibbs[i, i0:i1] = gibbs_values

            # Fill missing gibbs values (NaN entries) at start and end with extrapolated values
            g_row = gibbs[i]
            idx, = np.diff(np.isnan(g_row)).nonzero()
            if np.isnan(gibbs[i, 0]):
                fill_value = g_row[idx[0] + 1]
                g_row[:idx[0] + 1] = fill_value
            if np.isnan(g_row[-1]):
                fill_value = g_row[idx[-1]]
                g_row[idx[-1] + 1:] = fill_value

            X[i, 0: Npi, i0:i1] = kf.series.coverage.X
            timepoints[i, -Nti:] = kf.series.timepoints
            D[i, 0: Npi, -Nti:] = kf.series.uptake_corrected.T

        # Create pytorch tensors from input data, assign final shapes for matrix batch multiplication by tf.matmul
        dtype = torch.float64
        temperature_T = torch.tensor(temperature, dtype=dtype).unsqueeze(-1).unsqueeze(-1)  # Ns x 1 x 1
        k_int_T = torch.tensor(k_int, dtype=dtype).unsqueeze(-1)
        deltaG_T = torch.tensor(gibbs, dtype=dtype).unsqueeze(-1)
        timepoints_T = torch.tensor(timepoints, dtype=dtype).unsqueeze(1)
        X_T = torch.tensor(X, dtype=dtype)
        D_T = torch.tensor(D, dtype=dtype)

        deltaG_par = torch.nn.Parameter(deltaG_T)
        inputs = [temperature_T, X_T, k_int_T, timepoints_T]
        output_data = D_T

        #todo return as dict?
        return deltaG_par, inputs, output_data

    async def global_fit_async(self,  **kwargs):
        """see global fit"""
        client = await Client(self.cluster, asynchronous=True)
        fit_func = partial(self.global_fit, **kwargs)
        future = client.submit(fit_func)
        result = await future
        await client.close()

        return result

    def global_fit(self, r1=2, r2=5, epochs=100000, patience=50, stop_loss=0.05,
                   optimizer='SGD', **optimizer_kwargs):

        """

        Parameters
        ----------
        r1
        r2
        epochs
        patience
        stop_loss
        optimizer
        optimizer_kwargs

        Returns
        -------



        """
        # todo rewrite:
        # tensor_dict = self.setup_fit()
        # optimizer_kwargs, optimizer_klass =

        deltaG_par, inputs, output_data = self.setup_fit()
        model = DeltaGFit(deltaG_par)

        # Take default optimizer kwargs and update them with supplied kwargs
        optimizer_kwargs = {**optimizer_defaults.get(optimizer, {}), **optimizer_kwargs}

        optimizer_klass = getattr(torch.optim, optimizer)
        criterion = torch.nn.MSELoss(reduction='sum')

        #as property?
        def regularizer(param):
            d_ax1 = torch.abs(param[:, :-1, :] - param[:, 1:, :])
            d_ax2 = torch.abs(param - torch.mean(param, axis=0))
            reg_loss = r1 * torch.mean(d_ax1) + r2 * torch.mean(d_ax2)
            return reg_loss

        mse_loss, total_loss, returned_model = run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs,
                                                             model, criterion, regularizer, epochs=epochs,
                                                             patience=patience, stop_loss=stop_loss)

        result = TorchBatchFitResult(self, model, mse_loss=mse_loss, total_loss=total_loss)
        return result

    def global_fit_aligned(self, alignment_array, r1=2, r2=5, epochs=100000, patience=50, stop_loss=0.05,
                   optimizer='SGD', **optimizer_kwargs):
        #todo use run_optimizer function
        #todo (allow/force) indices as input

        """

        Parameters
        ----------
        alignment_array:
            Array of
        r1
        r2
        epochs
        patience
        stop_loss
        optimizer
        optimizer_kwargs

        Returns
        -------

        """
        r_numbers = np.cumsum(alignment_array != '-', axis=1)  # Residue numbers in alignment array
        aligned_bool = np.all(alignment_array != '-', axis=0)  # Array True where residues align
        aligned_residues = np.array([row[aligned_bool] for row in r_numbers])  # Residue numbers of aligned residues

        try:
            r_i0 = np.where(aligned_residues.min(axis=0) < self.interval[0])[0].max() + 1
        except ValueError:
            r_i0 = 0
        try:
            r_i1 = np.where(aligned_residues.max(axis=0) > self.interval[1] - 1)[0].min()
        except ValueError:
            r_i1 = aligned_residues.shape[1]

        # Crop aligned_residues to the part which corresponds to residues covered by measurements
        aligned_residues = aligned_residues[:, r_i0:r_i1]
        # Tranform residue number to index of corresponding deltaG values
        indices = aligned_residues - self.interval[0]

        i0 = torch.tensor(indices[0], dtype=torch.long)
        i1 = torch.tensor(indices[1], dtype=torch.long)
        deltaG_par, inputs, output_data = self.setup_fit()

        model = DeltaGFit(deltaG_par)

        #todo base class global fit function
        kwargs = optimizer_defaults.get(optimizer, {})
        kwargs.update(**optimizer_kwargs)

        optimizer_klass = getattr(torch.optim, optimizer)
        optimizer_obj = optimizer_klass(model.parameters(), **kwargs)

        criterion = torch.nn.MSELoss(reduction='sum')

        mse_loss = [torch.tensor(np.inf)]  # Mean squared loss only
        reg_loss = [torch.tensor(np.inf)]  # Loss including regularization loss
        stop = 0

        for epoch in range(epochs):
            optimizer_obj.zero_grad()
            output = model(*inputs)

            loss = criterion(output, output_data)
            mse_loss.append(loss)

            for pname, param in model.named_parameters():
                d_ax1 = torch.abs(param[:, :-1, :] - param[:, 1:, :])
                d_ax2 = torch.abs(param[0][i0] - param[1][i1])

                loss = loss + r1 * torch.mean(d_ax1) + r2 * torch.mean(d_ax2)

            reg_loss.append(loss)
            diff = reg_loss[-2] - loss
            if diff < stop_loss:
                stop += 1
                if stop > patience:
                    break
            else:
                stop = 0

            loss.backward()
            optimizer_obj.step()

        mse_loss = np.array([val.detach().numpy() for val in mse_loss])
        reg_loss = np.array([val.detach().numpy() for val in reg_loss])

        result = TorchBatchFitResult(self, model, mse_loss=mse_loss, total_loss=reg_loss)
        return result

    @property
    def temperature(self):
        return np.array([kf.temperature for kf in self.states])

    @property
    def exchanges(self):
        exchanges = np.zeros((self.Ns, self.Nr), dtype=bool)
        for i, kf in enumerate(self.states):
            interval_sample = kf.series.coverage.interval
            # Indices of residues
            i0 = interval_sample[0] - self.interval[0]
            i1 = interval_sample[1] - self.interval[0]
            exchanges[i, i0:i1] = kf.series.coverage['exchanges']

        return exchanges

    def do_guesses(self):
        raise NotImplementedError("This function should do guesses in batch on all kfs")




