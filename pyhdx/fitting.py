from scipy.optimize import fsolve
import numpy as np
from symfit import Fit, Variable, Parameter, exp, Model
from symfit.core.minimizers import DifferentialEvolution, MINPACK, LBFGSB, BFGS

r = Parameter('r', value=0.5)
tau1 = Parameter('tau1')
tau2 = Parameter('tau2')
t = Variable('t')
y = Variable('y')
model = Model({y: 100 *(1 - (r*exp(-t/tau1) + (1-r)*exp(-t/tau2)))})


def func_short(tau, tt, A):
    return 100 * (1 - np.exp(-tt / tau)) - A


def func_long(tau, tt, A, tau1):
    return 100 * (1 - (np.exp(-tt / tau1) + np.exp(-tt / tau))) - A


def initial_guess(t, d, r=0.5):
    tau1 = fsolve(func_short, 2, args=(t[2], d[2]))[0]
    tau2 = fsolve(func_long, 20, args=(t[-1], d[-1], tau1))[0]

    return tau1, tau2


def do_fitting(p_dict, n, times):
    scores_2d = np.stack([p_dict[n + '_' + str(t)].scores_average for t in times])
    s = p_dict[n + '_' + str[times[0]]]
    i = 0
    output = []
    for (j, state) in zip(s.cs, s.states):

        print(j)
        arr = scores_2d[:, i:j]
        i = j

        if np.any(np.isnan(arr)):  # states!
            output.append(np.nan)
            continue
        assert np.all(np.std(arr, axis=1) < 1e-10)

        d = arr[:, 0]
        chisq = np.inf
        # scan r for lowest chisquared
        for r_ in np.arange(0.1, 1, 0.1):
            t1, t2 = initial_guess(times, d, r=r_)
            f = model(np.array(times), r=r_, tau1=t1, tau2=t2)
            chisq_p = np.sum((d - f) ** 2)
            if chisq_p < chisq:
                t1f, t2f = t1, t2
                rf = r_
                chisq = chisq_p

        t1, t2 = initial_guess(times, d)

        tau1.value = t1
        tau2.value = min(t2, 100)
        r.value = r_

        fit = Fit(model, times, d, minimizer=MINPACK)
        res = fit.execute()

        rp = res.params['r'] * res.params['tau1'] + (1 - res.params['r']) * res.params['tau2']

        if np.isnan(rp):
            fit = Fit(model, times, d, minimizer=DifferentialEvolution)
            res = fit.execute()

            rp = res.params['r'] * res.params['tau1'] + (1 - res.params['r']) * res.params['tau2']
            if np.isnan(rp):
                print(j, state)
                print(d)
                raise ValueError('is still nan')
            
        output.append(rp)
