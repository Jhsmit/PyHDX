from pyhdx.support import get_reduced_blocks, temporary_seed
from pyhdx.models import Protein, HDXMeasurementSet
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
from dask.distributed import Client, worker_client
import dask
import warnings
import pandas as pd
from itertools import repeat


EmptyResult = namedtuple('EmptyResult', ['chi_squared', 'params'])
er = EmptyResult(np.nan, {k: np.nan for k in ['tau1', 'tau2', 'r']})


# ------------------------------------- #
# Rates fitting
# ------------------------------------- #

def get_bounds(times):
    """
    estimate default bound for rate fitting from a series of timepoints

    Parameters
    ----------
    times

    Returns
    -------

    """
    # todo document

    nonzero_times = times[np.nonzero(times)]
    t_first = np.min(nonzero_times)
    t_last = np.max(nonzero_times)
    b_upper = 50 * np.log(2) / t_first
    b_lower = np.log(2) / (t_last * 50)

    return b_lower, b_upper


def _prepare_wt_avg_fit(data_obj, model_type='association', bounds=None):
    series = data_obj  # todo refactor series
    bounds = bounds or get_bounds(data_obj.timepoints)

    arr = series.scores_stack.T  # data array
    i = 0
    # because intervals are inclusive, exclusive we need to add an extra entry to r_number for the final exclusive bound
    r_excl = np.append(series.coverage.r_number, [series.coverage.r_number[-1] + 1])

    models = []
    intervals = []  # Intervals; (start, end); (inclusive, exclusive)
    d_list = []
    for bl in series.coverage.block_length:
        d = arr[i]
        if np.all(np.isnan(d)):  # Skip non-coverage blocks
            i += bl
            continue
        intervals.append((r_excl[i], r_excl[i + bl]))
        d_list.append(d)
        if model_type == 'association':
            model = TwoComponentAssociationModel(bounds)
        elif model_type == 'dissociation':
            model = TwoComponentDissociationModel(bounds)
        else:
            raise ValueError('Invalid model type {}'.format(model_type))
        # res = EmptyResult(np.nan, {p.name: np.nan for p in model.sf_model.params})

        models.append(model)
        i += bl  # increment in block length does not equal move to the next start position

    return d_list, intervals, models


def fit_rates_half_time_interpolate(data_obj):
    """
    Calculates exchange rates based on weighted averaging followed by interpolation to determine half-time, which is
    then calculated to rates.

    Returns
    -------

    output: :class:`~np.ndarray`
        array with fields r_number, rate

    """
    # todo this is uing the soon to be depcrecated coverage object
    interpolated = np.array(
        [np.interp(50, d_uptake, data_obj.timepoints) for d_uptake in data_obj.scores_stack.T])

    output = np.empty_like(interpolated, dtype=[('r_number', int), ('rate', float)])
    output['r_number'] = data_obj.coverage.r_number
    output['rate'] = np.log(2) / interpolated

    protein = Protein(output, index='r_number')
    t50FitResult = namedtuple('t50FitResult', ['output'])

    result = t50FitResult(output=protein)

    return result


def fit_rates_weighted_average(data_obj, bounds=None, chisq_thd=20, model_type='association', client=None, pbar=None):
    """
    Block length _should_ be equal to the block length of all measurements in the series, provided that all coverage
    is the same

    Parameters
    ----------
    chisq_max

    Returns
    -------

    """

    d_list, intervals, models = _prepare_wt_avg_fit(data_obj, model_type=model_type, bounds=bounds)
    if pbar:
        raise NotImplementedError()
    else:
        inc = lambda: None

    results = []

    if client is None:
        for d, model in zip(d_list, models):
            result = fit_kinetics(data_obj.timepoints, d, model, chisq_thd=chisq_thd)
            results.append(result)
    else:
        iterables = [[data_obj.timepoints]*len(d_list), d_list, models]

        if isinstance(client, Client):
            futures = client.map(fit_kinetics, *iterables, chisq_thd=chisq_thd)
            results = client.gather(futures)
        elif client == 'worker_client':
            with worker_client() as client:
                futures = client.map(fit_kinetics, *iterables, chisq_thd=chisq_thd)
                results = client.gather(futures)


    fit_result = KineticsFitResult(data_obj, intervals, results, models)

    return fit_result


def fit_rates(data_obj, method='wt_avg', **kwargs):
    """
    Fit observed rates of exchange to HDX-MS data in `data_obj`

    Parameters
    ----------
    data_obj: KineticsSeries
    method: :obj:`str`
        Method to use to determine rates of exchange
    kwargs
        Additional kwargs passed to fitting

    Returns
    -------

    fit_result : class;KinetcisFitresult

    """

    if method == 'wt_avg':
        result = fit_rates_weighted_average(data_obj, **kwargs)
    else:
        raise ValueError(f"Invalid value for 'method': {method}")

    return result


def fit_kinetics(t, d, model, chisq_thd=100):
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

#
# def fit_global(data, model):
#     fit = Fit(model.sf_model, **data)
#     res = fit.execute()
#     return res


# ------------------------------------- #
# Gibbs fitting
# ------------------------------------- #

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


def regularizer_1d(r1, param):
    return r1 * torch.mean(torch.abs(param[:-1] - param[1:]))


def regularizer_2d(r1, r2, param):
    #todo allow regularization wrt reference rather than mean
    d_ax1 = torch.abs(param[:, :-1, :] - param[:, 1:, :])
    d_ax2 = torch.abs(param - torch.mean(param, axis=0))
    reg_loss = r1 * torch.mean(d_ax1) + r2 * torch.mean(d_ax2)
    return reg_loss


def regularizer_2d_aligned(r1, r2, indices, param):
    i0 = indices[0]
    i1 = indices[1]
    d_ax1 = torch.abs(param[:, :-1, :] - param[:, 1:, :])
    d_ax2 = torch.abs(param[0][i0] - param[1][i1])

    reg_loss = r1 * torch.mean(d_ax1) + r2 * torch.mean(d_ax2)
    return reg_loss


def fit_gibbs_global(data_object, initial_guess, r1=2, epochs=100000, patience=50, stop_loss=0.05,
               optimizer='SGD', **optimizer_kwargs):
    #todo @tejas: Missing docstring
    """Pytorch global fitting"""

    tensors = data_object.get_tensors()
    inputs = [tensors[key] for key in ['temperature', 'X', 'k_int', 'timepoints']]
    output_data = tensors['uptake']

    if isinstance(initial_guess, pd.Series):
        initial_guess = initial_guess.to_numpy()

    assert len(initial_guess) == data_object.Nr, "Invalid length of initial guesses"
    #todo dtype config
    dtype = torch.float64
    deltaG_par = torch.nn.Parameter(torch.tensor(initial_guess, dtype=dtype).unsqueeze(-1))  #reshape (nr, 1)
    #deltaG_par = torch.nn.Parameter(torch.Tensor(initial_guess).unsqueeze(-1))

    model = DeltaGFit(deltaG_par)
    criterion = torch.nn.MSELoss(reduction='sum')

    # Take default optimizer kwargs and update them with supplied kwargs
    optimizer_kwargs = {**optimizer_defaults.get(optimizer, {}), **optimizer_kwargs}  # Take defaults and override with user-specified
    optimizer_klass = getattr(torch.optim, optimizer)

    reg_func = partial(regularizer_1d, r1)

    # returned_model is the same object as model
    mse_loss, total_loss, returned_model = run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs,
                                                         model, criterion, reg_func, epochs=epochs,
                                                         patience=patience, stop_loss=stop_loss)

    result = TorchSingleFitResult(data_object, model,
                                  mse_loss=mse_loss, total_loss=total_loss)

    return result


def fit_gibbs_global_batch(hdx_set, initial_guess, r1=2, r2=5, epochs=100000, patience=50, stop_loss=0.05,
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
    # todo still some repeated code with fit_gibbs single
    tensors = hdx_set.get_tensors()
    inputs = [tensors[key] for key in ['temperature', 'X', 'k_int', 'timepoints']]
    output_data = tensors['uptake']

    assert initial_guess.shape == (hdx_set.Ns, hdx_set.Nr), "Invalid shape of initial guesses"

    dtype = torch.float64
    deltaG_par = torch.nn.Parameter(torch.tensor(initial_guess, dtype=dtype).reshape(hdx_set.Ns, hdx_set.Nr, 1))

    model = DeltaGFit(deltaG_par)
    criterion = torch.nn.MSELoss(reduction='sum')

    # Take default optimizer kwargs and update them with supplied kwargs
    optimizer_kwargs = {**optimizer_defaults.get(optimizer, {}), **optimizer_kwargs}  # Take defaults and override with user-specified
    optimizer_klass = getattr(torch.optim, optimizer)

    reg_func = partial(regularizer_2d, r1, r2)
    mse_loss, total_loss, returned_model = run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs,
                                                         model, criterion, reg_func, epochs=epochs,
                                                         patience=patience, stop_loss=stop_loss)

    result = TorchBatchFitResult(hdx_set, model, mse_loss=mse_loss, total_loss=total_loss)
    return result


def fit_gibbs_global_batch_aligned(hdx_set, initial_guess, r1=2, r2=5, epochs=100000, patience=50, stop_loss=0.05,
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

    assert hdx_set.Ns == 2, 'Aligned batch fitting is limited to two states'

    tensors = hdx_set.get_tensors()
    inputs = [tensors[key] for key in ['temperature', 'X', 'k_int', 'timepoints']]
    output_data = tensors['uptake']

    assert initial_guess.shape == (hdx_set.Ns, hdx_set.Nr), "Invalid shape of initial guesses"

    dtype = torch.float64
    deltaG_par = torch.nn.Parameter(torch.tensor(initial_guess, dtype=dtype).reshape(hdx_set.Ns, hdx_set.Nr, 1))

    model = DeltaGFit(deltaG_par)
    criterion = torch.nn.MSELoss(reduction='sum')

    # Take default optimizer kwargs and update them with supplied kwargs
    optimizer_kwargs = {**optimizer_defaults.get(optimizer, {}), **optimizer_kwargs}  # Take defaults and override with user-specified
    optimizer_klass = getattr(torch.optim, optimizer)

    if hdx_set.aligned_indices is None:
        raise ValueError("No alignment added to HDX measurements")

    indices = [torch.tensor(i, dtype=torch.long) for i in hdx_set.aligned_indices]

    reg_func = partial(regularizer_2d_aligned, r1, r2, indices)
    mse_loss, total_loss, returned_model = run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs,
                                                         model, criterion, reg_func, epochs=epochs,
                                                         patience=patience, stop_loss=stop_loss)

    result = TorchBatchFitResult(hdx_set, model, mse_loss=mse_loss, total_loss=total_loss)
    return result


"""
this might still serve some use
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
"""


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

