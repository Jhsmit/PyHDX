from pyhdx import Coverage

from scipy.optimize import fsolve
import numpy as np
from symfit import Fit, Variable, Parameter, exp, Model
from symfit.core.minimizers import DifferentialEvolution, Powell
from collections import namedtuple


#module level (non dummy) parameters are likely to cause all sorts of problems
r = Parameter('r', value=0.5, min=0, max=1)
tau1 = Parameter('tau1', min=0, max=5)
tau2 = Parameter('tau2', min=0, max=100)
t = Variable('t')
y = Variable('y')
model = Model({y: 100 *(1 - (r*exp(-t/tau1) + (1-r)*exp(-t/tau2)))})


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


def initial_guess(t, d):
    """
    Calculates initial guesses for fitting of two-component kinetic uptake reaction

    Parameters
    ----------
    t : :class:~`numpy.ndarray`
        Array with time points
    d : :class:~`numpy.ndarray`
        Array with uptake values

    Returns
    -------
    tau1 : :obj:`float`
        Initial guess for short time component
    tau2 : :obj:`float`
        Initial guess for long time component

    """
    tau1 = fsolve(func_short, 2, args=(t[2], d[2]))[0]
    tau2 = fsolve(func_long, 20, args=(t[-2], d[-2], tau1))[0]

    return tau1, tau2


EmptyResult = namedtuple('EmptyResult', ['chi_squared', 'params'])
er = EmptyResult(np.nan, {k: np.nan for k in ['tau1', 'tau2', 'r']})


def fit_kinetics(t, d, chisq_thd):
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
        return er

    t1, t2 = initial_guess(t, d)

    tau1.value = t1
    tau2.value = min(t2, 200)
    r.value = 0.5

    fit = Fit(model, t, d, minimizer=Powell)
    res = fit.execute()
    rp = res.params['r'] * res.params['tau1'] + (1 - res.params['r']) * res.params['tau2']

    if np.isnan(rp) or res.chi_squared > chisq_thd or res.params['r'] > 1 or res.params['r'] < 0:

        #TODO add thread lock here
        fit = Fit(model, t, d, minimizer=DifferentialEvolution)
        res = fit.execute(workers=-1)

    return res


class KineticsFitting(object):

    def __init__(self, k_series):
        #todo check if coverage
        self.k_series = k_series


      #  self.coverage = Coverage(k_series[0].data)


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
        i = 0
        for j, d in enumerate(arr):
            i += 1
            if j == len(arr) - 1:  # End of array
                block_length.append(i)
                continue
            elif np.all(d == arr[j + 1]):  # Next element is the same
                continue
            else:
                block_length.append(i)
                i = 0

                res = fit_kinetics(self.k_series.times, d, chisq_thd=chisq_thd)
                results.append(res)

        return results, block_length





