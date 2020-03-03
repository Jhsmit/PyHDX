from pyhdx import Coverage

from scipy.optimize import fsolve
import numpy as np
from symfit import Fit, Variable, Parameter, exp, Model
from symfit.core.minimizers import DifferentialEvolution, Powell
from collections import namedtuple

from tqdm.auto import tqdm
#
# #module level (non dummy) parameters are likely to cause all sorts of problems
# r = Parameter('r', value=0.5, min=0, max=1)
# tau1 = Parameter('tau1', min=0, max=5)
# tau2 = Parameter('tau2', min=0, max=100)
# t = Variable('t')
# y = Variable('y')
# model = Model({y: 100 *(1 - (r*exp(-t/tau1) + (1-r)*exp(-t/tau2)))})
#


class KineticsModel(object):
    par_index = 0
    var_index = 0

    def __init__(self):
        self.names = {}
        self.sf_model = None

    def make_parameter(self, name, value=1., min=None, max=None):
        dummy_name = 'pyhdx_par_{}'.format(self.par_index)
        KineticsModel.par_index += 1
        p = Parameter(dummy_name, value=value, min=min, max=max)
        self.names[name] = dummy_name
        return p

    def make_variable(self, name):
        dummy_name = 'pyhdx_var_{}'.format(self.var_index)
        KineticsModel.var_index += 1
        v = Variable(dummy_name)
        self.names[name] = dummy_name
        return v

    @property
    def r_names(self):
        """ dictionary of dummy_name: name"""
        return {v: k for k, v in self.names.items()}

    def get_parameter(self, name):
        """returns the parameter object with human-readable name name"""
        dummy_name = self.names[name]
        par_names = list([p.name for p in self.sf_model.params])
        idx = par_names.index(dummy_name)
        return self.sf_model.params[idx]


class BiexpIncreaseModel(KineticsModel):
    """Two Phase Association"""
    def __init__(self):
        super(BiexpIncreaseModel, self).__init__()

        r = self.make_parameter('r', value=0.5, min=0, max=1)
        tau1 = self.make_parameter('tau1', min=0, max=5)
        tau2 = self.make_parameter('tau2', min=0, max=100)
        t = self.make_variable('t')
        y = self.make_variable('y')

        self.sf_model = Model({y: 100 * (1 - (r * exp(-t / tau1) + (1 - r) * exp(-t / tau2)))})

    def initial_guess(self, t, d):
        """
        Calculates initial guesses for fitting of two-component kinetic uptake reaction

        Parameters
        ----------
        t : :class:~`numpy.ndarray`
            Array with time points
        d : :class:~`numpy.ndarray`
            Array with uptake values

        """
        tau1_v = fsolve(func_short, 2, args=(t[2], d[2]))[0]
        tau2_v = fsolve(func_long, 20, args=(t[-2], d[-2], tau1_v))[0]

        tau1_p = self.get_parameter('tau1')
        tau1_p.value = tau1_v
        tau2_p = self.get_parameter('tau2')
        tau2_p.value = tau2_v
        r_p = self.get_parameter('r')
        r_p.value = 0.5

    def get_tau(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """

        r = params[self.names['r']]
        tau1 = params[self.names['tau1']]
        tau2 = params[self.names['tau2']]

        tau = r * tau1 + (1 - r) * tau2

        return tau

    def get_rate(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """
        tau = self.get_tau(**params)
        return 1/tau


class MonoexpIncreaseModel(KineticsModel):
    """One Phase Association"""
    def __init__(self):
        super(MonoexpIncreaseModel, self).__init__()
        tau1 = self.make_parameter('tau1', min=0, max=100)
        t = self.make_variable('t')
        y = self.make_variable('y')

        self.sf_model = Model({y: 100 * (1 - exp(-t / tau1))})

    def initial_guess(self, t, d):
        """
        Calculates initial guesses for fitting of two-component kinetic uptake reaction

        Parameters
        ----------
        t : :class:~`numpy.ndarray`
            Array with time points
        d : :class:~`numpy.ndarray`
            Array with uptake values

        """
        tau1_v = fsolve(func_short, 2, args=(t[3], d[3]))[0]

        tau1_p = self.get_parameter('tau1')
        tau1_p.value = tau1_v

    def get_tau(self, **params):
        tau1 = params[self.names['tau1']]
        return tau1

    def get_rate(self, **params):
        tau = self.get_tau(**params)
        return 1/tau


def func_short(tau, tt, A):
    """
    Function to estimate the short time component

    Parameters
    ----------
    tau : :obj:`float`
        Lifetime
    tt : :obj:`float`
        Selected time point
    A : :obj:`float`
        Target amplitude

    Returns
    -------
    A_t : :obj:`float`
        Amplitude difference given tau, tt, A

    """
    return 100 * (1 - np.exp(-tt / tau)) - A


def func_long(tau, tt, A, tau1):
    """
    Function to estimate the short time component

    Parameters
    ----------
    tau : :obj:`float`
        Lifetime
    tt : :obj:`float`
        Selected time point
    A : :obj:`float`
        Target amplitude
    tau1: : obj:`float`
        Lifetime of short time component

    Returns
    -------
    A_t : :obj:`float`
        Amplitude difference given tau, tt, A, tau1

    """
    return 100 * (1 - (0.5 * np.exp(-tt / tau1) + 0.5 * np.exp(-tt / tau))) - A


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
    chisq_thd: :obj:`float`
        Threshold chi squared above which the fitting is repeated with the Differential Evolution algorithm.

    Returns
    -------
    res : :class:`~symfit.FitResults`
        Symfit fitresults object.
    """
    if np.any(np.isnan(d)):  # states!
        er = EmptyResult(np.nan, {p.name: np.nan for p in model.sf_model.params})
        return er

    model.initial_guess(t, d)
    fit = Fit(model.sf_model, t, d, minimizer=Powell)
    res = fit.execute()
    rate = model.get_rate(**res.params)
    try:
        r = res.params[model.names['r']]
    except KeyError:
        r = 1

    if np.isnan(rate) or res.chi_squared > chisq_thd or r > 1 or r < 0:
        print(res.chi_squared)
        #TODO add thread lock here
        fit = Fit(model.sf_model, t, d, minimizer=DifferentialEvolution)
        res = fit.execute(workers=-1)

    return res


class KineticsFitting(object):

    def __init__(self, k_series):
        #todo check if coverage is the same!
        self.k_series = k_series
        self.result = None

    @property
    def scores_stack(self):
        """uptake scores to fit in a 2d stack"""
        scores_2d = np.stack([v.scores_average for v in self.k_series])
        return scores_2d

    @property
    def scores_norm(self):
    # Normalized to 100 array of scores
        scores_norm = 100 * (self.scores_stack / self.scores_stack[-1, :][np.newaxis, :])
        return scores_norm

    @property
    def scores_peptides(self):
        scores_peptides = np.stack([v.scores for v in self.k_series])
        return scores_peptides


    def global_fitting(self):
        """
        fit (per section) in time

        Returns
        -------

        """

    def do_fitting(self, chisq_thd=20):
        """
        Block length _should_ be equal to the block length of all measurements in the series, provided that all coverage
        is the same

        Parameters
        ----------
        chisq_max

        Returns
        -------

        """

        arr = self.scores_norm.T
        block_length = []
        results = []
        models = []
        i = 0
        for j, d in enumerate(tqdm(arr)):
            i += 1
            # End of array, or new block approaching,
            if j == len(arr) - 1 or not np.allclose(d, arr[j + 1], equal_nan=True):
                block_length.append(i)
                i = 0

                model = BiexpIncreaseModel()
                res = fit_kinetics(self.k_series.times, d, model, chisq_thd=chisq_thd)
                results.append(res)
                models.append(model)

        self.result = KineticsFitResult(results, models, block_length)
        return self.result


class KineticsFitResult(object):
    """
    this fit results is only for wt avg fitting
    """
    def __init__(self, results, models, block_length):
        assert len(results) == len(models)
        assert len(models) == len(block_length)

        self.results = results
        self.models = models
        self.block_length = block_length

    def __len__(self):
        return len(self.results)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.results[item], self.models[item], self.block_length[item]
        else:
            return KineticsFitResult(self.results[item], self.models[item], self.block_length[item])

    def __iter__(self):
        iterable = [(r, m, b) for r, m, b in zip(self.results, self.models, self.block_length)]
        return iterable.__iter__()

    def get_param(self, name, repeat=True):
        """
        Get an array of parameter with name `name` from the fit result. The length of the array is equal to the
        number of amino acids.

        Parameters
        ----------
        name : :obj:`str`
            Name of the parameter to extract
        repeat : :obj:`bool`
            If true the parameter array is repeated by block size

        Returns
        -------
        par_arr : :class:`~numpy.ndarray`
            Array with parameter values

        """

        names = [model.names[name] for model in self.models]
        arr = np.array([res.params[name] for res, name in zip(self.results, names)])
        if repeat:
            par_arr = np.repeat(arr, self.block_length)
        else:
            par_arr = arr
        return par_arr

    @property
    def rate(self):
        """Returns an array with the exchange rates"""
        rates = np.array([model.get_rate(**res.params) for model, res in zip(self.models, self.results)])
        return np.repeat(rates, self.block_length)




