from pyhdx import Coverage

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

    def fine_fitting_fixed_blocks(self, initial_result):
        """ initial_result: KineticsFitResult object from global_fitting"""

        assert self.k_series.uniform
        split = self.k_series.split()
        print(len(split))

        total_cs = np.append(self.k_series.cov.start, self.k_series.cov.start + np.cumsum(self.k_series.cov.block_length))
        rate_output = np.empty(len(self.k_series.cov.r_number))
        rate_output[:] = np.nan
        for section in split.values():  # Section block is 7 one block
            s, e = section.cov.start, section.cov.end
            i0, i1 = np.searchsorted(total_cs, [s, e + 1])
            # section_result = initial_result[i0:i1]

            #on fail, this section has block length 8 while section result has 1, 7

            kf = KineticsFitting(section)  #deeeeep. dont do fine_fitting on this guy
            print('make model')
            #model = LSQKineticsDoubleConstantBlocks(initial_result, kf)
            model = LSQKineticsDoubleCustomBlocks(initial_result, kf)


            print('make data dict')
            data_dict = {d_var.name: scores for d_var, scores in zip(model.d_vars, kf.scores_peptides.T)}
            data_dict[model.t_var.name] = section.times
            print('make fit')
            fit = Fit(model.sf_model, **data_dict)
            print('do fit')
            result = fit.execute()
            print('done')
            rate = model.get_rate(**result.params)
            i0, i1 = np.searchsorted(self.k_series.cov.r_number, [s, e])
            rate_output[i0:i1+1] = rate

        return self.k_series.cov.r_number, rate_output

    def fine_fitting(self, initial_result):
        """ initial_result: KineticsFitResult object from global_fitting"""

        assert self.k_series.uniform
        split = self.k_series.split()
        print(len(split))

        total_cs = np.append(self.k_series.cov.start, self.k_series.cov.start + np.cumsum(self.k_series.cov.block_length))
        rate_output = np.empty(len(self.k_series.cov.r_number))
        rate_output[:] = np.nan
        for section in split.values():  # Section block is 7 one block
            s, e = section.cov.start, section.cov.end
            i0, i1 = np.searchsorted(total_cs, [s, e + 1])
            section_result = initial_result[i0:i1]

            #on fail, this section has block length 8 while section result has 1, 7

            kf = KineticsFitting(section)  #deeeeep. dont do fine_fitting on this guy
            print('make model')
            model = LSQKinetics(section_result, kf)

            print('make data dict')
            data_dict = {d_var.name: scores for d_var, scores in zip(model.d_vars, kf.scores_peptides.T)}
            data_dict[model.t_var.name] = section.times
            print('make fit')
            fit = Fit(model.sf_model, **data_dict)
            print('do fit')
            result = fit.execute()
            print('done')
            rate = model.get_rate(**result.params)
            i0, i1 = np.searchsorted(self.k_series.cov.r_number, [s, e])
            rate_output[i0:i1+1] = rate

        return self.k_series.cov.r_number, rate_output

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


class LSQKineticsDoubleCustomBlocks(KineticsModel): #TODO find a better name (lstsq)
    #todo block length is redundant
    def __init__(self, initial_result, kf_section):  #todo does this really need to be kf but can be series?
        """

        Parameters
        ----------
        initial_result array with r_number and rate
        kf_section kineticsfitting object for the section
        """
        super(LSQKineticsDoubleCustomBlocks, self).__init__()
        # print(block_lengths)
        # total = np.sum(block_lengths)
        t_var = self.make_variable('t')

        #Assemble terms spanning accross blocks of residues
        #blocks with only 1 residue have only 1 time component
        terms = []
        r_number = initial_result['r_number']

        block_length = list(kf_section.k_series.cov.block_length.copy())


        # Merge blocks
        max_size = 2
        i = 0
        subblock = 0
        while i < len(block_length):
            curr_size = block_length[i]
            if curr_size <= max_size:
                del block_length[i]
                subblock += curr_size
            elif subblock != 0:
                # if block_length > 1:
                block_length.insert(i, subblock)
                subblock = 0
                i += 1
            else:
                i += 1
        if subblock != 0:
            block_length.append(subblock)

        blocks = block_length
        print(blocks)

        changes = True
        min_size = 5
        while changes:
            i = 0
            changes = False
            while i < len(block_length):
                curr_size = block_length[i]
                if curr_size < min_size:
                    changes = True
                    if i == 0:  # beginning of list
                        block_length[i + 1] += curr_size
                        del block_length[i]
                    elif i == len(block_length) - 1:  # end of the list
                        block_length[i - 1] += curr_size
                        del block_length[i]
                    else:
                        if block_length[i - 1] < block_length[i + 1]:
                            block_length[i - 1] += curr_size
                        else:
                            block_length[i + 1] += curr_size
                        del block_length[i]
                else:
                    i += 1

        print(block_length)
        r_index = kf_section.k_series.cov.start
        for i, bl in enumerate(blocks):
            current = r_index + (bl // 2)
            idx = np.searchsorted(r_number, current)
            t1v = 1 / initial_result['rate1'][idx]
            t2v = 1 / initial_result['rate2'][idx]
            r_index += bl

            # t1v = value / 5
            # t2v = value * 5
            rv = 0.5
            tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1/80, value=t1v)
            tau2 = self.make_parameter('tau2_{}'.format(i), max=30, min=1/80, value=t2v)
            r = self.make_parameter('r_{}'.format(i), max=1, min=0, value=rv)
            term = (1 - (r * exp(-t_var / tau1) + (1 - r) * exp(-t_var / tau2)))
            terms.append(term)

        #Iterate over rows (peptides) and add terms together which make one peptide
        model_dict = {}
        d_vars = []
        cs = np.insert(np.cumsum(blocks) + kf_section.k_series.cov.start, 0, kf_section.k_series.cov.start)
        for i, entry in enumerate(kf_section.k_series.cov.data):
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

        self.block_length = blocks

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


class LSQKineticsDoubleConstantBlocks(KineticsModel): #TODO find a better name (lstsq)
    #todo block length is redundant
    def __init__(self, initial_result, kf_section):  #todo does this really need to be kf but can be series?
        """

        Parameters
        ----------
        initial_result array with r_number and rate
        kf_section kineticsfitting object for the section
        """
        super(LSQKineticsDoubleConstantBlocks, self).__init__()
        # print(block_lengths)
        # total = np.sum(block_lengths)
        t_var = self.make_variable('t')

        #Assemble terms spanning accross blocks of residues
        #blocks with only 1 residue have only 1 time component
        terms = []

        block_size = 10
        initial_block = 5
        num_repeats = (kf_section.k_series.cov.prot_len - initial_block) // block_size
        remainder = (kf_section.k_series.cov.prot_len - initial_block) % block_size

        r_number = initial_result['r_number']
        rate = initial_result['rate']

        blocks = [initial_block] + [block_size] * num_repeats
        if remainder:
            blocks += [remainder]
        r_index = kf_section.k_series.cov.start
        for i, bl in enumerate(blocks):
            current = r_index + (bl // 2)
            idx = np.searchsorted(r_number, current)
            t1v = 1 / initial_result['rate1'][idx]
            t2v = 1 / initial_result['rate2'][idx]
            r_index += bl

            # t1v = value / 5
            # t2v = value * 5
            rv = 0.5
            tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1/40, value=t1v)
            tau2 = self.make_parameter('tau2_{}'.format(i), max=30, min=1/40, value=t2v)
            r = self.make_parameter('r_{}'.format(i), max=1, min=0, value=rv)
            term = (1 - (r * exp(-t_var / tau1) + (1 - r) * exp(-t_var / tau2)))
            terms.append(term)

        #Iterate over rows (peptides) and add terms together which make one peptide
        model_dict = {}
        d_vars = []
        cs = np.insert(np.cumsum(blocks) + kf_section.k_series.cov.start, 0, kf_section.k_series.cov.start)
        for i, entry in enumerate(kf_section.k_series.cov.data):
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

        self.block_length = blocks

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


class LSQKineticsConstantBlocks(KineticsModel): #TODO find a better name (lstsq)
    #todo block length is redundant
    def __init__(self, initial_result, kf_section):  #todo does this really need to be kf but can be series?
        """

        Parameters
        ----------
        initial_result array with r_number and rate
        kf_section kineticsfitting object for the section
        """
        super(LSQKineticsConstantBlocks, self).__init__()
        # print(block_lengths)
        # total = np.sum(block_lengths)
        t_var = self.make_variable('t')

        #Assemble terms spanning accross blocks of residues
        #blocks with only 1 residue have only 1 time component
        terms = []

        block_size = 5
        initial_block = 2
        num_repeats = (kf_section.k_series.cov.prot_len - initial_block) // block_size
        remainder = (kf_section.k_series.cov.prot_len - initial_block) % block_size

        r_number = initial_result['r_number']
        rate = initial_result['rate']

        blocks = [initial_block] + [block_size] * num_repeats
        if remainder:
            blocks += [remainder]
        r = kf_section.k_series.cov.start
        for i, bl in enumerate(blocks):
            current = r + (bl // 2)
            idx = np.searchsorted(r_number, current)
            value = 1 / rate[idx]
            r += bl
            print('val', value)
            tau1 = self.make_parameter('tau_{}'.format(i), max=30, min=1 / 40, value=value)
            term = (1 - exp(-t_var / tau1))
            terms.append(term)

        #Iterate over rows (peptides) and add terms together which make one peptide
        model_dict = {}
        d_vars = []
        cs = np.insert(np.cumsum(blocks) + kf_section.k_series.cov.start, 0, kf_section.k_series.cov.start)
        for i, entry in enumerate(kf_section.k_series.cov.data):
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

        self.block_length = blocks

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
            tau = params[self.names['tau_{}'.format(i)]]
            tau_list += [tau]*bl

        print(tau_list)
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



class LSQKinetics(KineticsModel): #TODO find a better name (lstsq)
    #todo block length is redundant
    def __init__(self, initial_result, kf_section):
        """

        Parameters
        ----------
        initial_result kineticsresult object from initial fitting
        kf_section kineticsfitting object for the section
        """
        super(LSQKinetics, self).__init__()
        # print(block_lengths)
        # total = np.sum(block_lengths)
        t_var = self.make_variable('t')

        #Assemble terms spanning accross blocks of residues
        #blocks with only 1 residue have only 1 time component
        terms = []
        for i, (r, m, bl) in enumerate(initial_result):
            if bl == 1:
                t1v = r.params[m.names['tau1']]
                t2v = r.params[m.names['tau2']]
                rv = r.params[m.names['r']]
                value = rv * t1v + (1 - rv) * t2v
                tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1 / 40, value=value)

                term = (1 - exp(-t_var / tau1))
            else:
                t1v = r.params[m.names['tau1']]
                t2v = r.params[m.names['tau2']]
                rv = r.params[m.names['r']]
                tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1 / 40, value=t1v)
                tau2 = self.make_parameter('tau2_{}'.format(i), max=30, min=1 / 40, value=t2v)
                r = self.make_parameter('r_{}'.format(i), max=1, min=0, value=rv)

                term = (1 - (r*exp(-t_var / tau1) + (1-r)*exp(-t_var/tau2)))
            terms.append(term)

        #Iterate over rows (peptides) and add terms together which make one peptide
        model_dict = {}
        d_vars = []
        for i, x_row in enumerate(kf_section.k_series.cov.X_red_norm):
            d_var = self.make_variable('d_{}'.format(i))
            d_vars.append(d_var)
            rhs = reduce(add, [100*fraction*term for fraction, term in zip(x_row, terms)])
            model_dict[d_var] = rhs

        self.d_vars = d_vars
        self.t_var = t_var
        self.sf_model = CallableModel(model_dict)

        self.block_length = initial_result.block_length

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
            if bl == 1:
                tau1 = params[self.names['tau1_{}'.format(i)]]
                tau_list.append(tau1)
            else:
                tau1 = params[self.names['tau1_{}'.format(i)]]
                tau2 = params[self.names['tau2_{}'.format(i)]]
                r = params[self.names['r_{}'.format(i)]]

                tau = r * tau1 + (1 - r) * tau2
                tau_list += [tau]*bl

        print(tau_list)
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
