from collections import namedtuple
from functools import partial

import numpy as np
import pandas as pd
import torch
from dask.distributed import Client, worker_client
from symfit import Fit
from symfit.core.minimizers import DifferentialEvolution, Powell
from tqdm import trange

from pyhdx.fit_models import SingleKineticModel, TwoComponentAssociationModel, TwoComponentDissociationModel
from pyhdx.fitting_torch import DeltaGFit, TorchSingleFitResult, TorchBatchFitResult, TORCH_DTYPE, TORCH_DEVICE
from pyhdx.models import Protein
from pyhdx.support import temporary_seed

EmptyResult = namedtuple('EmptyResult', ['chi_squared', 'params'])
er = EmptyResult(np.nan, {k: np.nan for k in ['tau1', 'tau2', 'r']})


# Reguarlizers act on ΔG values, which are in kJ/mol and range typically from 0 to 40000 J/mol.
# Therefore they are scaled by a factor 10000 such that they are near one (as D values are also near one)
REGULARIZATION_SCALING = 1e-4


# default values

PATIENCE = 50
STOP_LOSS = 1e-5
EPOCHS = 100000
R1 = 1
R2 = 1

optimizer_defaults = {
    'SGD': {
        'lr': 1e4,
        'momentum': 0.5,
        'nesterov': True
    },
}

# ------------------------------------- #
# Rates fitting
# ------------------------------------- #

def get_bounds(times):
    """
    estimate default bound for rate fitting from a series of timepoints

    Parameters
    ----------
    times : array_like

    Returns
    -------
    bounds : :obj:`tuple`
        lower and upper bounds

    """
    # todo document

    nonzero_times = times[np.nonzero(times)]
    t_first = np.min(nonzero_times)
    t_last = np.max(nonzero_times)
    b_upper = 50 * np.log(2) / t_first
    b_lower = np.log(2) / (t_last * 50)

    return b_lower, b_upper


def _prepare_wt_avg_fit(hdxm, model_type='association', bounds=None):
    """

    Parameters
    ----------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
    model_type
    bounds tuple

    Returns
    -------



    """
    bounds = bounds or get_bounds(hdxm.timepoints)

    arr = hdxm.rfu_residues  # Data array
    i = 0
    # because intervals are inclusive, exclusive we need to add an extra entry to r_number for the final exclusive bound
    r_excl = np.append(hdxm.coverage.r_number, [hdxm.coverage.r_number[-1] + 1])

    models = []
    intervals = []  # Intervals; (start, end); (inclusive, exclusive)
    d_list = []
    for bl in hdxm.coverage.block_length:
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


def fit_rates_half_time_interpolate(hdxm):
    """
    Calculates exchange rates based on weighted averaging followed by interpolation to determine half-time, which is
    then calculated to rates.

    Parameters
    ----------

    hdxm : :class:`~pyhdx.models.HDXMeasurement`


    Returns
    -------

    output: dataclass
        dataclass with fit result

    """
    interpolated = np.array(
        [np.interp(0.5, d_uptake, hdxm.timepoints) for d_uptake in hdxm.rfu_residues])

    output = np.empty_like(interpolated, dtype=[('r_number', int), ('rate', float)])
    output['r_number'] = hdxm.coverage.r_number
    output['rate'] = np.log(2) / interpolated

    protein = Protein(output, index='r_number')
    t50FitResult = namedtuple('t50FitResult', ['output']) # todo dataclass?

    result = t50FitResult(output=protein)

    return result


def fit_rates_weighted_average(hdxm, bounds=None, chisq_thd=0.20, model_type='association', client=None, pbar=None):
    """
    Fit a model specified by 'model_type' to D-uptake kinetics. D-uptake is weighted averaged across peptides per
    timepoint to obtain residue-level D-uptake.

    Parameters
    ----------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
    bounds : :obj:`tuple`, optional
        Tuple of lower and upper bounds of rate constants in the model used.
    chisq_thd : :obj:`float`
        Threshold of chi squared result, values above will trigger a second round of fitting using DifferentialEvolution
    model_type : :obj:`str`
        Missing docstring
    client : : ??
        Controls delegation of fitting tasks to Dask clusters. Options are: `None`: Do not use task, fitting is done
        in the local thread in a for loop. :class: Dask Client : Uses the supplied Dask client to schedule fitting task.
        `worker_client`: The function was ran by a Dask worker and the additional fitting tasks created are scheduled
        on the same Cluster.
    pbar:
        Not implemented

    Returns
    -------

    fit_result : :class:`~pyhdx.fitting.KineticsFitResult`

    """
    d_list, intervals, models = _prepare_wt_avg_fit(hdxm, model_type=model_type, bounds=bounds)
    if pbar:
        raise NotImplementedError()
    else:
        inc = lambda: None

    results = []

    if client is None:
        for d, model in zip(d_list, models):
            result = fit_kinetics(hdxm.timepoints, d, model, chisq_thd=chisq_thd)
            results.append(result)
    else:
        iterables = [[hdxm.timepoints]*len(d_list), d_list, models]

        if isinstance(client, Client):
            futures = client.map(fit_kinetics, *iterables, chisq_thd=chisq_thd)
            results = client.gather(futures)
        elif client == 'worker_client':
            with worker_client() as client:
                futures = client.map(fit_kinetics, *iterables, chisq_thd=chisq_thd)
                results = client.gather(futures)

    fit_result = KineticsFitResult(hdxm, intervals, results, models)

    return fit_result


def fit_rates(hdxm, method='wt_avg', **kwargs):
    """
    Fit observed rates of exchange to HDX-MS data in `hdxm`

    Parameters
    ----------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
    method : :obj:`str`
        Method to use to determine rates of exchange
    kwargs
        Additional kwargs passed to fitting

    Returns
    -------

    fit_result : :class:`~pyhdx.fitting.KineticsFitResult`

    """

    if method == 'wt_avg':
        result = fit_rates_weighted_average(hdxm, **kwargs)
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
    model : :class:`~pyhdx.fit_models.KineticsModel`
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
    """ Check if the obtained fit result is within bounds"""
    for param in fit_result.model.params:
        value = fit_result.params[param.name]
        if value < param.min:
            return False
        elif value > param.max:
            return False
    return True


# ------------------------------------- #
# Gibbs fitting
# ------------------------------------- #


def run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs, model, criterion, regularizer,
                  epochs=EPOCHS, patience=PATIENCE, stop_loss=STOP_LOSS, callbacks=None, tqdm=True):
    """

    Runs optimization/fitting of PyTorch model.

    Parameters
    ----------
    inputs : :obj:`list`
        List of input Tensors
    output_data : :class:`~torch.Tensor`
        comparison data to model output
    optimizer_klass : :mod:`~torch.optim`
    optimizer_kwargs : :obj:`dict`
        kwargs to pass to pytorch optimizer
    model : :class:`~torch.nn.Module`
        pytorch model
    criterion: callable
        loss function
    regularizer callable
        regularizer function
    epochs : :obj:`int`
        Max number of epochs
    patience : :obj:`int`
        Number of epochs with less progress than `stop_loss` before terminating optimization
    stop_loss : :obj:`float`
        Threshold of optimization value below which no progress is made
    callbacks: :obj:`list` or `None`
        List of callback functions
    tqdm : :obj:`bool`
        Toggle tqdm progress bar

    Returns
    -------

    """

    optimizer_obj = optimizer_klass(model.parameters(), **optimizer_kwargs)

    #todo these seeds should be temporary
    np.random.seed(43)
    torch.manual_seed(43)

    callbacks = callbacks or []
    losses_list = [[np.inf]]

    def closure():
        output = model(*inputs)
        loss = criterion(output, output_data)
        losses_list.append([loss.item()])  # store mse loss
        reg_loss_tuple = regularizer(model.deltaG)
        for r in reg_loss_tuple:
            loss += r

        losses_list[-1] += [r.item() for r in reg_loss_tuple]  # store reg losses

        loss.backward()
        return loss

    stop = 0
    iter = trange(epochs) if tqdm else range(epochs)
    for epoch in iter:
        optimizer_obj.zero_grad()
        loss = optimizer_obj.step(closure)

        for cb in callbacks:
            cb(epoch, model, optimizer_obj)

        diff = sum(losses_list[-2]) - sum(losses_list[-1])
        if diff < stop_loss:
            stop += 1
            if stop > patience:
                break
        else:
            stop = 0

    return np.array(losses_list[1:]), model


def regularizer_1d(r1, param):
    reg_loss = r1 * torch.mean(torch.abs(param[:-1] - param[1:]))
    return reg_loss * REGULARIZATION_SCALING,


def regularizer_2d_mean(r1, r2, param):
    #todo allow regularization wrt reference rather than mean
    #param shape: Ns x Nr x 1
    d_ax1 = torch.abs(param[:, :-1, :] - param[:, 1:, :])
    d_ax2 = torch.abs(param - torch.mean(param, axis=0))

    return r1 * torch.mean(d_ax1) * REGULARIZATION_SCALING, r2 * torch.mean(d_ax2) * REGULARIZATION_SCALING


def regularizer_2d_reference(r1, r2, param):
    d_ax1 = torch.abs(param[:, :-1, :] - param[:, 1:, :])
    d_ax2 = torch.abs(param - param[0])[1:]

    return r1 * torch.mean(d_ax1) * REGULARIZATION_SCALING, r2 * torch.mean(d_ax2) * REGULARIZATION_SCALING


def regularizer_2d_aligned(r1, r2, indices, param):
    i0 = indices[0]
    i1 = indices[1]
    d_ax1 = torch.abs(param[:, :-1, :] - param[:, 1:, :])
    d_ax2 = torch.abs(param[0][i0] - param[1][i1])

    return r1 * torch.mean(d_ax1) * REGULARIZATION_SCALING, r2 * torch.mean(d_ax2) * REGULARIZATION_SCALING


def _loss_df(losses_array):
    """transforms losses array to losses dataframe
    first column in losses array is mse loss, rest are regularzation losses
    """

    loss_df = pd.DataFrame(losses_array, columns=['mse_loss'] + [f'reg_{i + 1}' for i in range(losses_array.shape[1] - 1)])
    loss_df.index.name = 'epoch'
    loss_df.index += 1

    return loss_df


def fit_gibbs_global(hdxm, initial_guess, r1=R1, epochs=EPOCHS, patience=PATIENCE, stop_loss=STOP_LOSS,
                     optimizer='SGD', callbacks=None, **optimizer_kwargs):
    """
    Fit Gibbs free energies globally to all D-uptake data in the supplied hdxm

    Parameters
    ----------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
    initial_guess : :class:`~pandas.Series` or :class:`~numpy.ndarray`
        Gibbs free energy initial guesses (shape Nr)
    r1 : :obj:`float`
    epochs
    patience
    stop_loss
    optimizer : :obj:`str`
    callbacks: :obj:`list` or None
        List of callback objects. call signature is cb(epoch, model, optimizer)
    optimizer_kwargs

    Returns
    -------

    """
    #todo @tejas: Missing docstring
    """Pytorch global fitting"""

    fit_keys = ['r1', 'epochs', 'patience', 'stop_loss', 'optimizer']
    locals_dict = locals()
    fit_kwargs = {k: locals_dict[k] for k in fit_keys}

    tensors = hdxm.get_tensors()
    inputs = [tensors[key] for key in ['temperature', 'X', 'k_int', 'timepoints']]
    output_data = tensors['uptake']

    if isinstance(initial_guess, pd.Series):
        initial_guess = initial_guess.to_numpy()

    assert len(initial_guess) == hdxm.Nr, "Invalid length of initial guesses"

    dtype = torch.float64
    deltaG_par = torch.nn.Parameter(torch.tensor(initial_guess, dtype=TORCH_DTYPE, device=TORCH_DEVICE).unsqueeze(-1))  #reshape (nr, 1)

    model = DeltaGFit(deltaG_par)
    criterion = torch.nn.MSELoss(reduction='mean')

    # Take default optimizer kwargs and update them with supplied kwargs
    optimizer_kwargs = {**optimizer_defaults.get(optimizer, {}), **optimizer_kwargs}  # Take defaults and override with user-specified
    optimizer_klass = getattr(torch.optim, optimizer)

    reg_func = partial(regularizer_1d, r1)

    # returned_model is the same object as model
    losses_array, returned_model = run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs,
                                                         model, criterion, reg_func, epochs=epochs,
                                                         patience=patience, stop_loss=stop_loss, callbacks=callbacks)
    losses = _loss_df(losses_array)
    fit_kwargs.update(optimizer_kwargs)
    result = TorchSingleFitResult(hdxm, model, losses=losses, **fit_kwargs)

    return result


def fit_gibbs_global_batch(hdx_set, initial_guess, r1=R1, r2=R2, r2_reference=False, epochs=EPOCHS, patience=PATIENCE, stop_loss=STOP_LOSS,
               optimizer='SGD', callbacks=None, **optimizer_kwargs):
    """
    Batch fit gibbs free energies to multiple HDX measurements

    Parameters
    ----------
    hdx_set : :class:`~pyhdx.models.HDXMeasurementSet`
    initial_guess
    r1
    r2
    r2_reference,
        if true the first dataset is used as a reference to calculate r2 differences, otherwise the mean is used
    epochs
    patience
    stop_loss
    optimizer
    callbacks: :obj:`list` or None
        List of callback objects. call signature is cb(epoch, model, optimizer)
    optimizer_kwargs

    Returns
    -------

    """
    # todo still some repeated code with fit_gibbs single

    fit_keys = ['r1', 'r2', 'r2_reference', 'epochs', 'patience', 'stop_loss', 'optimizer', 'callbacks']
    locals_dict = locals()
    fit_kwargs = {k: locals_dict[k] for k in fit_keys}

    if r2_reference:
        reg_func = partial(regularizer_2d_reference, r1, r2)
    else:
        reg_func = partial(regularizer_2d_mean, r1, r2)

    return _batch_fit(hdx_set, initial_guess, reg_func, fit_kwargs, optimizer_kwargs)


def fit_gibbs_global_batch_aligned(hdx_set, initial_guess, r1=R1, r2=R2, epochs=EPOCHS, patience=PATIENCE,
                                   stop_loss=STOP_LOSS, optimizer='SGD', callbacks=None, **optimizer_kwargs):
    """
    Batch fit gibbs free energies to two HDX measurements. The supplied HDXMeasurementSet must have alignment information
    (supplied by HDXMeasurementSet.add_alignment)


    Parameters
    ----------
    hdx_set : :class:`~pyhdx.models.HDXMeasurementSet`
    initial_guess
    r1
    r2
    epochs
    patience
    stop_loss
    callbacks: :obj:`list` or None
        List of callback objects. call signature is cb(epoch, model, optimizer)
    optimizer
    optimizer_kwargs

    Returns
    -------

    """

    assert hdx_set.Ns == 2, 'Aligned batch fitting is limited to two states'
    if hdx_set.aligned_indices is None:
        raise ValueError("No alignment added to HDX measurements")

    indices = [torch.tensor(i, dtype=torch.long) for i in hdx_set.aligned_indices]
    reg_func = partial(regularizer_2d_aligned, r1, r2, indices)

    fit_keys = ['r1', 'r2', 'epochs', 'patience', 'stop_loss', 'optimizer', 'callbacks']
    locals_dict = locals()
    fit_kwargs = {k: locals_dict[k] for k in fit_keys}

    return _batch_fit(hdx_set, initial_guess, reg_func, fit_kwargs, optimizer_kwargs)


def _batch_fit(hdx_set, initial_guess, reg_func, fit_kwargs, optimizer_kwargs):
    # @tejas docstrings
    tensors = hdx_set.get_tensors()
    inputs = [tensors[key] for key in ['temperature', 'X', 'k_int', 'timepoints']]
    output_data = tensors['uptake']

    assert initial_guess.shape == (hdx_set.Ns, hdx_set.Nr), "Invalid shape of initial guesses"

    deltaG_par = torch.nn.Parameter(torch.tensor(initial_guess, dtype=TORCH_DTYPE, device=TORCH_DEVICE).reshape(hdx_set.Ns, hdx_set.Nr, 1))

    model = DeltaGFit(deltaG_par)
    criterion = torch.nn.MSELoss(reduction='mean')

    # Take default optimizer kwargs and update them with supplied kwargs
    optimizer_kwargs = {**optimizer_defaults.get(fit_kwargs['optimizer'], {}), **optimizer_kwargs}  # Take defaults and override with user-specified
    optimizer_klass = getattr(torch.optim, fit_kwargs['optimizer'])

    loop_kwargs = {k: fit_kwargs[k] for k in ['epochs', 'patience', 'stop_loss']}
    loop_kwargs['callbacks'] = fit_kwargs.pop('callbacks')
    losses_array, returned_model = run_optimizer(inputs, output_data, optimizer_klass, optimizer_kwargs,
                                                         model, criterion, reg_func, **loop_kwargs)
    losses = _loss_df(losses_array)
    fit_kwargs.update(optimizer_kwargs)
    result = TorchBatchFitResult(hdx_set, model, losses=losses, **fit_kwargs)

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
    Fit result object. Generally used for initial guess results.

    Parameters
    ----------

    hdxm : :class:`~pyhdx.models.HDXMeasurement`
    intervals:
    results:
    models:

    """
    def __init__(self, hdxm, intervals, results, models):
        """
        each model with corresponding interval covers a region in the protein corresponding to r_number
        """
        assert len(results) == len(models)
#        assert len(models) == len(block_length)
        self.hdxm = hdxm
        self.r_number = hdxm.coverage.r_number
        self.intervals = intervals  #inclusive, excluive
        self.results = results
        self.models = models

    @property
    def model_type(self):
        # Most likely all current instances have model_type `single`
        if np.all([isinstance(m, SingleKineticModel) for m in self.models]):
            return 'Single'
        else:
            raise ValueError('Unsupported model types')

    def __call__(self, timepoints):
        """call the result with timepoints to get fitted uptake per peptide back"""
        #todo outdated
        d_list = []
        if self.model_type == 'Single':
            for time in timepoints:
                p = self.get_p(time)
                p = np.nan_to_num(p)
                d = self.hdxm.coverage.X.dot(p)
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
        #todo pd series/dataframe in favour of searchsorted
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

