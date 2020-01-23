from scipy.optimize import fsolve
import numpy as np
from symfit import Fit, Variable, Parameter, exp, Model
from symfit.core.minimizers import DifferentialEvolution, Powell
from collections import namedtuple


#module level parameters are likely to cause all sort of problems
r = Parameter('r', value=0.5, min=0, max=1)
tau1 = Parameter('tau1', min=0, max=5)
tau2 = Parameter('tau2', min=0, max=100)
t = Variable('t')
y = Variable('y')
model = Model({y: 100 *(1 - (r*exp(-t/tau1) + (1-r)*exp(-t/tau2)))})


def func_short(tau, tt, A):
    return 100 * (1 - np.exp(-tt / tau)) - A


def func_long(tau, tt, A, tau1):
    return 100 * (1 - (0.5 * np.exp(-tt / tau1) + 0.5 * np.exp(-tt / tau))) - A


def initial_guess(t, d):
    tau1 = fsolve(func_short, 2, args=(t[2], d[2]))[0]
    tau2 = fsolve(func_long, 20, args=(t[-2], d[-2], tau1))[0]

    return tau1, tau2


EmptyResult = namedtuple('EmptyResult', ['chi_squared', 'params'])
er = EmptyResult(np.nan, {k: np.nan for k in ['tau1', 'tau2', 'r']})


def fit_kinetics(t, d, chisq_thd):
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
        fit = Fit(model, t, d, minimizer=DifferentialEvolution)
        res = fit.execute(workers=-1)

    return res
