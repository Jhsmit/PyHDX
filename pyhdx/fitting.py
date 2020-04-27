from pyhdx import Coverage
from pyhdx.support import get_constant_blocks, get_reduced_blocks

from scipy.optimize import fsolve
import numpy as np
from symfit import Fit, Variable, Parameter, exp, Model, CallableModel
from symfit.core.minimizers import DifferentialEvolution, Powell
from collections import namedtuple
from functools import reduce
from operator import add

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
        self.names = {}  # human name: dummy name
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


class TwoComponentAssociationModel(KineticsModel):
    """Two componenent Association"""
    def __init__(self):
        super(TwoComponentAssociationModel, self).__init__()

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

    def initial_grid(self, t, d, step=15):
        kmax = 5 * np.log(1-0.98) / -t[1]
        d_final = np.min([0.95, d[-1]/100])  # todo refactor norm
        kmin = np.log(1-d_final) / -t[-1]

        tau_space = np.logspace(np.log10(1/kmax), np.log10(1/kmin), num=step, endpoint=True)
        r_space = np.linspace(0.05, 0.95, num=step, endpoint=True)

        guess = np.column_stack([tau_space, tau_space, r_space])
        return guess

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


class OneComponentAssociationModel(KineticsModel):
    """One component Association"""
    def __init__(self):
        super(OneComponentAssociationModel, self).__init__()
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
        raise ValueError('There shouldnt be NaNs anymore')
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
        #TODO add thread lock here
        fit = Fit(model.sf_model, t, d, minimizer=DifferentialEvolution)
        #grid = model.initial_grid(t, d, step=5)
        res = fit.execute()

    return res


class KineticsFitting(object):

    def __init__(self, k_series):
        #todo check if coverage is the same!
        self.k_series = k_series
        self.result = None

    @property
    def scores_stack(self):
        raise DeprecationWarning()
        print('move to series')
        # todo move this to series
        """uptake scores to fit in a 2d stack"""
        scores_2d = np.stack([v.scores_average for v in self.k_series])
        return scores_2d

    @property
    def scores_norm(self):
    # Normalized to 100 array of scores
        raise DeprecationWarning()
        print('where is this used?')
        scores_norm = 100 * (self.scores_stack / self.scores_stack[-1, :][np.newaxis, :])
        return scores_norm

    @property
    def scores_peptides(self):
        raise DeprecationWarning()
        print('move to series')
        scores_peptides = np.stack([v.scores for v in self.k_series])
        return scores_peptides

    def lsq_fit_blocks(self, initial_result, model='bi', block_func=get_reduced_blocks, **block_kwargs):
        """ initial_result: KineticsFitResult object from global_fitting"""

        assert self.k_series.uniform
        split = self.k_series.split()

        rate_output = np.empty(len(self.k_series.cov.r_number))
        rate_output[:] = np.nan

        results = []
        models = []
        intervals = []
        for section in split.values():  # Section block is 7 one block
            s, e = section.cov.start, section.cov.end  #inclusive, inclusive (source file)
            intervals.append((s, e + 1))  #inclusive, exclusive

            blocks = block_func(section, **block_kwargs)
            model = LSQKinetics(initial_result, section, blocks, model=model)

            # scores_peptides = np.stack([v.scores for v in section])
            # sp1 = section.scores_peptides
            # assert np.all(scores_peptides == sp1)
            data_dict = {d_var.name: scores for d_var, scores in zip(model.d_vars, section.scores_peptides.T)}
            data_dict[model.t_var.name] = section.times
            fit = Fit(model.sf_model, **data_dict)
            result = fit.execute()

            results.append(result)
            models.append(model)

        fit_result = KineticsFitResult(self.k_series.cov.r_number, intervals, results, models)

        return fit_result

    def weighted_avg_fit(self, chisq_thd=20):
        """
        Block length _should_ be equal to the block length of all measurements in the series, provided that all coverage
        is the same

        Parameters
        ----------
        chisq_max

        Returns
        -------

        """

        results = []
        models = []
        intervals = [] # Intervals; (start, end); (inclusive, exclusive)

        for k, series in self.k_series.split().items():
            arr = series.scores_stack.T
            #arr = self.scores_stack.T  # todo use nonnormed scores as they should be normalized already
            i = 0
            for bl in series.cov.block_length:
                intervals.append((series.cov.start + i, series.cov.start + i + bl + 1))
                d = arr[i]
                model = TwoComponentAssociationModel()
                res = fit_kinetics(series.times, d, model, chisq_thd=chisq_thd)
                results.append(res)
                models.append(model)
                i += bl

        self.result = KineticsFitResult(self.k_series.cov.r_number, intervals, results, models)
        return self.result


class KineticsFitResult(object):
    """
    this fit results is only for wt avg fitting
    """
    def __init__(self, r_number, intervals, results, models):
        assert len(results) == len(models)
#        assert len(models) == len(block_length)
        self.r_number = r_number
        self.intervals = intervals
        self.results = results
        self.models = models


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
    def tau(self):
        """Returns an array with the exchange rates"""
        output = np.full_like(self.r_number, np.nan, dtype=float)
        for (s, e), result, model in zip(self.intervals, self.results, self.models):
            rate = model.get_tau(**result.params)
            i0, i1 = np.searchsorted(self.r_number, [s, e])
            output[i0:i1] = rate
        return output

    @property
    def rate(self):
        """Returns an array with the exchange rates"""
        return 1 / self.tau

        # output = np.full_like(self.r_number, np.nan, dtype=float)
        # for (s, e), result, model in zip(self.intervals, self.results, self.models):
        #     rate = model.get_rate(**result.params)
        #     i0, i1 = np.searchsorted(self.r_number, [s, e])
        #     output[i0:i1] = rate
        # return output

    def get_output(self, names):
        dtype = [('r_number', int)] + [(name, float) for name in names]
        array = np.full_like(self.r_number, np.nan, dtype=dtype)
        array['r_number'] = self.r_number
        for name in names:
            try:
                array[name] = getattr(self, name)
            except AttributeError:
                array[name] = self.get_param(name)
        return array


class LSQKinetics(KineticsModel): #TODO find a better name (lstsq)
    #todo block length is redundant
    #what is keeping it as attribute?
    #not really because its now needed to keep as we're not storing it anymore on FitREsults
    def __init__(self, initial_result, k_series, blocks, model='bi'):
    #todo allow for setting of min/max value of parameters
        """

        Parameters
        ----------
        initial_result array with r_number and rate
        k_series kineticsseries object for the section
        """
        super(LSQKinetics, self).__init__()
        t_var = self.make_variable('t')

        self.block_length = blocks
        self.model = model
        terms = []

        r_number = initial_result['r_number']

        r_index = k_series.cov.start
        for i, bl in enumerate(blocks):
            current = r_index + (bl // 2)
            idx = np.searchsorted(r_number, current)
            if model == 'mono' or (model == 'mixed' and bl == 1):
                print('mono')
                #TODO update names
                value = 1 / initial_result['rate'][idx]  #this should be the average of the block range
                r_index += bl
                tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1 / 70, value=value)
                term = (1 - exp(-t_var / tau1))
                terms.append(term)
            else:
                # t1v = initial_result['tau1'][idx]  #TODO update names
                # t2v = initial_result['tau2'][idx]
                r_index += bl


                t1v = min(initial_result['tau1'][idx], initial_result['tau2'][idx])
                t2v = max(initial_result['tau1'][idx], initial_result['tau2'][idx])
                print(t1v, t2v)


                rv = 0.5
                tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1 / 70, value=t1v)
                tau2 = self.make_parameter('tau2_{}'.format(i), max=30, min=1 / 70, value=t2v)
                r = self.make_parameter('r_{}'.format(i), max=1, min=0, value=rv)
                term = (1 - (r * exp(-t_var / tau1) + (1 - r) * exp(-t_var / tau2)))
                terms.append(term)

        #Iterate over rows (peptides) and add terms together which make one peptide
        model_dict = {}
        d_vars = []
        cs = np.insert(np.cumsum(blocks) + k_series.cov.start, 0, k_series.cov.start)
        for i, entry in enumerate(k_series.cov.data):
            d_var = self.make_variable('d_{}'.format(i))
            d_vars.append(d_var)

            s, e = entry['start'], entry['end']
            length = e - s + 1
            used_blocks = np.diff(np.clip((cs - s).copy(), a_min=0, a_max=length))
            fractions = used_blocks / length

            rhs = reduce(add, [100*fraction*term for fraction, term in zip(fractions, terms)])
            model_dict[d_var] = rhs

        self.d_vars = d_vars
        self.t_var = t_var
        self.sf_model = CallableModel(model_dict)

    def get_param_values(self, name, **params):
        """returns a list of parameters with name name which should have been indexed parameters
        params repeat during blocks

        """

        values = [params[self.names[f'{name}_{i}']] for i, _ in enumerate(self.block_length)]
        return np.repeat(values, self.block_length)

    def get_tau(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """

        tau_list = []
        for i, bl in enumerate(self.block_length):
            if self.model == 'mono' or (self.model == 'mixed' and bl == 1):
                tau = params[self.names['tau1_{}'.format(i)]]
            else:
                tau1 = params[self.names['tau1_{}'.format(i)]]
                tau2 = params[self.names['tau2_{}'.format(i)]]
                r = params[self.names['r_{}'.format(i)]]

                tau = r * tau1 + (1 - r) * tau2
            tau_list += [tau] * bl

        return np.array(tau_list)

    def get_tau(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """

        tau_list = []
        for i, bl in enumerate(self.block_length):
            if self.model == 'mono' or (self.model == 'mixed' and bl == 1):
                tau = params[self.names['tau_{}'.format(i)]]
            else:
                tau1 = params[self.names['tau1_{}'.format(i)]]
                tau2 = params[self.names['tau2_{}'.format(i)]]
                r = params[self.names['r_{}'.format(i)]]

                tau = r * tau1 + (1 - r) * tau2
            tau_list += [tau] * bl

        return np.array(tau_list)

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


