import numpy as np
from symfit import Parameter, Variable, Model, exp
from scipy.optimize import fsolve


class KineticsModel(object):
    """
    Base class for kinetics models. Main function is to generate :ref:`symfit` Variables and Parameters. The class
    attributes `par_index` and `var_index` are used to make sure names used by :ref:`symfit` are unique and their
    mapping to user-defined names are stored in the `names` dictionary.

    Parameters
    ----------

    bounds : :obj:`tuple`
        Tuple of default `min`, `max` parameters to use.

    Attributes
    ----------

    names : :obj:`dict`
        Dictionary which maps human-readable names (keys) to dummy names (values)
    sf_model : :class:`~symfit.Model`
        The `symfit` model which describes this model. Implemented by subclasses.

    """

    par_index = 0
    var_index = 0

    def __init__(self, bounds):
        if bounds[1] < bounds[0]:
            raise ValueError('Lower bound must be smaller than upper bound')
        self.bounds = bounds
        self.names = {}  # human name: dummy name
        self.sf_model = None

    def make_parameter(self, name, value=None, min=None, max=None):
        """
        Create a new :class:`~symfit.Parameter`.

        Parameters
        ----------
        name : :obj:`str`
            Human-readable name for the parameter
        value : :obj:`float`
            Initial guess value
        min : :obj:`float`
            Lower bound value. If `None`, the value from `bounds` is used.
        max : :obj:`float`
            Lower bound value. If `None`, the value from `bounds` is used.

        Returns
        -------
        p : :class:`~symfit.Parameter`

        """
        min = min if min is not None else self.bounds[0]
        max = max if max is not None else self.bounds[1]

        value = value or np.mean(self.bounds)
        dummy_name = 'pyhdx_par_{}'.format(self.par_index)
        KineticsModel.par_index += 1
        p = Parameter(dummy_name, value=value, min=min, max=max)
        self.names[name] = dummy_name
        return p

    def make_variable(self, name):
        """
        Create a new :class:`~symfit.Variable`.

        Parameters
        ----------
        name : :obj:`str`
            Human-readable name for the variable

        Returns
        -------
        p : :class:`~symfit.Variable`

        """
        dummy_name = 'pyhdx_var_{}'.format(self.var_index)
        KineticsModel.var_index += 1
        v = Variable(dummy_name)
        self.names[name] = dummy_name
        return v

    @property
    def r_names(self):
        """:obj:`dict`: Reverse dictionary of the variable and parameter names"""
        return {v: k for k, v in self.names.items()}

    def get_parameter(self, name):
        """
        Get the parameter with the Human-readable name `name`

        Parameters
        ----------
        name : :obj:`str`
            Name of the parameter to retrieve

        Returns
        -------
        parameter : :class:`~symfit.Parameter`

        """
        dummy_name = self.names[name]
        par_names = list([p.name for p in self.sf_model.params])
        idx = par_names.index(dummy_name)
        parameter = self.sf_model.params[idx]
        return parameter


class SingleKineticModel(KineticsModel):
    """
    Base class for models which fit only a single set (slice) of time, uptake points
    """


class TwoComponentAssociationModel(SingleKineticModel):
    """Two componenent Association"""
    def __init__(self, bounds):
        super(TwoComponentAssociationModel, self).__init__(bounds)

        r = self.make_parameter('r', value=0.5, min=0, max=1)
        k1 = self.make_parameter('k1')
        k2 = self.make_parameter('k2')
        t = self.make_variable('t')
        y = self.make_variable('y')

        self.sf_model = Model({y: (1 - (r * exp(-k1*t) + (1 - r) * exp(-k2*t)))})

    def __call__(self, t, **params):
        """call model at time t, returns uptake values of peptides"""
        time_var = self.names['t']
        params[time_var] = t
        return self.sf_model(**params)

    def initial_guess(self, t, d):
        """
        Calculates initial guesses for fitting of two-component kinetic uptake reaction

        Parameters
        ----------
        t : :class:`~numpy.ndarray`
            Array with time points
        d : :class:`~numpy.ndarray`
            Array with uptake values

        """
        k1_v = fsolve(func_short_ass, 1 / 2, args=(t[2], d[2]))[0]
        k2_v = fsolve(func_long_ass, 1 / 20, args=(t[-2], d[-2], k1_v))[0]

        k1_p = self.get_parameter('k1')
        k1_p.value = k1_v
        k2_p = self.get_parameter('k2')
        k2_p.value = k2_v
        r_p = self.get_parameter('r')
        r_p.value = 0.5

    def initial_grid(self, t, d, step=15):
        kmax = 5 * np.log(1-0.98) / -t[1]
        d_final = np.min([0.95, d[-1]])  # todo refactor norm
        kmin = np.log(1-d_final) / -t[-1]

        tau_space = np.logspace(np.log10(1/kmax), np.log10(1/kmin), num=step, endpoint=True)
        r_space = np.linspace(0.05, 0.95, num=step, endpoint=True)

        guess = np.column_stack([tau_space, tau_space, r_space])
        return guess

    def get_rate(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """

        #todo generalize duplicate code for other models
        r = params[self.names['r']]
        k1 = params[self.names['k1']]
        k2 = params[self.names['k2']]

        #tau = r * tau1 + (1 - r) * tau2

        x0 = k1 if r > 0.5 else k2

        t_log = fsolve(self.min_func, np.log(x0), args=(k1, k2, r))
        try:
            assert np.round(self.min_func(t_log, k1, k2, r), 5) == 0, 'Failed to find half life root'
        except AssertionError:
            print('uhoh')
            print(k1, k2, r)
            print(np.round(self.min_func(t_log, k1, k2, r), 5) == 0)
        k = np.asscalar(np.log(2) / np.exp(t_log))

        return k

    def get_tau(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """
        k = self.get_rate(**params)
        return 1/k

    @staticmethod
    def min_func(t_log, k1, k2, r):
        t = np.exp(t_log)
        return 0.5 - r * np.exp(-k1*t) - (1 - r) * np.exp(-k2*t)


class OneComponentAssociationModel(SingleKineticModel):
    """One component Association"""
    def __init__(self, bounds):
        super(OneComponentAssociationModel, self).__init__(bounds)
        k1 = self.make_parameter('k1')
        t = self.make_variable('t')
        y = self.make_variable('y')

        self.sf_model = Model({y: (1 - exp(-k1*t))})

    def __call__(self, t, **params):
        """call model at time t, returns uptake values of peptides"""

        time_var = self.names['t']
        params[time_var] = t
        return self.sf_model(**params)

    def initial_guess(self, t, d):
        """
        Calculates initial guesses for fitting of two-component kinetic uptake reaction

        Parameters
        ----------
        t : :class:`~numpy.ndarray`
            Array with time points
        d : :class:`~numpy.ndarray`
            Array with uptake values

        """
        k1_v = fsolve(func_short_ass, 1 / 2, args=(t[3], d[3]))[0]

        k1_p = self.get_parameter('k1')
        k1_p.value = k1_v

    def get_rate(self, **params):
        k1 = params[self.names['k1']]
        return k1

    def get_tau(self, **params):
        k = self.get_rate(**params)
        return 1/k


class TwoComponentDissociationModel(SingleKineticModel):
    """Two componenent Association"""
    def __init__(self, bounds):
        super(TwoComponentDissociationModel, self).__init__(bounds)

        r = self.make_parameter('r', value=0.5, min=0, max=1)
        k1 = self.make_parameter('k1')
        k2 = self.make_parameter('k2')
        t = self.make_variable('t')
        y = self.make_variable('y')

        self.sf_model = Model({y: (r * exp(-k1*t) + (1 - r) * exp(-k2*t))})

    def __call__(self, t, **params):
        """call model at time t, returns uptake values of peptides"""
        time_var = self.names['t']
        params[time_var] = t
        return self.sf_model(**params)

    def initial_guess(self, t, d):
        """
        Calculates initial guesses for fitting of two-component kinetic uptake reaction

        Parameters
        ----------
        t : :class:`~numpy.ndarray`
            Array with time points
        d : :class:`~numpy.ndarray`
            Array with uptake values

        """
        k1_v = fsolve(func_short_ass, 1 / 2, args=(t[2], d[2]))[0]
        k2_v = fsolve(func_long_ass, 1 / 20, args=(t[-2], d[-2], k1_v))[0]

        k1_p = self.get_parameter('k1')
        k1_p.value = k1_v
        k2_p = self.get_parameter('k2')
        k2_p.value = k2_v
        r_p = self.get_parameter('r')
        r_p.value = 0.5

    def initial_grid(self, t, d, step=15):
        kmax = 5 * np.log(1-0.98) / -t[1]
        d_final = np.min([0.95, d[-1]])  # todo refactor norm
        kmin = np.log(1-d_final) / -t[-1]

        tau_space = np.logspace(np.log10(1/kmax), np.log10(1/kmin), num=step, endpoint=True)
        r_space = np.linspace(0.05, 0.95, num=step, endpoint=True)

        guess = np.column_stack([tau_space, tau_space, r_space])
        return guess

    def get_rate(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """

        #todo generalize duplicate code for other models
        r = params[self.names['r']]
        k1 = params[self.names['k1']]
        k2 = params[self.names['k2']]

        #tau = r * tau1 + (1 - r) * tau2

        x0 = k1 if r > 0.5 else k2

        t_log = fsolve(self.min_func, np.log(x0), args=(k1, k2, r))
        try:
            assert np.round(self.min_func(t_log, k1, k2, r), 5) == 0, 'Failed to find half life root'
        except AssertionError:
            print('uhoh')
            print(k1, k2, r)
            print(np.round(self.min_func(t_log, k1, k2, r), 5) == 0)
        k = np.asscalar(np.log(2) / np.exp(t_log))

        return k

    def get_tau(self, **params):
        """

        Parameters
        ----------
        params

        key value where keys are the dummy names

        Returns
        -------

        """
        k = self.get_rate(**params)
        return 1/k

    @staticmethod
    def min_func(t_log, k1, k2, r):
        t = np.exp(t_log)
        return - 0.5 + r * np.exp(-k1*t) + (1 - r) * np.exp(-k2*t)


class OneComponentDissociationModel(SingleKineticModel):
    """One component Association"""
    def __init__(self, bounds):
        super(OneComponentDissociationModel, self).__init__(bounds)
        k1 = self.make_parameter('k1')
        t = self.make_variable('t')
        y = self.make_variable('y')

        self.sf_model = Model({y: exp(-k1*t)})

    def __call__(self, t, **params):
        """call model at time t, returns uptake values of peptides"""

        time_var = self.names['t']
        params[time_var] = t
        return self.sf_model(**params)

    def initial_guess(self, t, d):
        """
        Calculates initial guesses for fitting of two-component kinetic uptake reaction

        Parameters
        ----------
        t : :class:`~numpy.ndarray`
            Array with time points
        d : :class:`~numpy.ndarray`
            Array with uptake values

        """
        k1_v = fsolve(func_short_ass, 1 / 2, args=(t[3], d[3]))[0]

        k1_p = self.get_parameter('k1')
        k1_p.value = k1_v

    def get_rate(self, **params):
        k1 = params[self.names['k1']]
        return k1

    def get_tau(self, **params):
        k = self.get_rate(**params)
        return 1/k


def func_short_dis(k, tt, A):
    """
    Function to estimate the fast time component

    Parameters
    ----------
    k : :obj:`float`
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
    return np.exp(-k * tt) - A


def func_long_dis(k, tt, A, k1):
    """
    Function to estimate the short time component

    Parameters
    ----------
    k : :obj:`float`
        rate
    tt : :obj:`float`
        Selected time point
    A : :obj:`float`
        Target amplitude
    k1 : :obj:`float`
        Rate of fast time component

    Returns
    -------
    A_t : :obj:`float`
        Amplitude difference given tau, tt, A, tau1

    """
    return (0.5 * np.exp(-k1*tt) + 0.5 * np.exp(-k*tt)) - A


def func_short_ass(k, tt, A):
    """
    Function to estimate the fast time component

    Parameters
    ----------
    k : :obj:`float`
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
    return (1 - np.exp(-k * tt)) - A


def func_long_ass(k, tt, A, k1):
    """
    Function to estimate the short time component

    Parameters
    ----------
    k : :obj:`float`
        rate
    tt : :obj:`float`
        Selected time point
    A : :obj:`float`
        Target amplitude
    k1 : :obj:`float`
        Rate of fast time component

    Returns
    -------
    A_t : :obj:`float`
        Amplitude difference given tau, tt, A, tau1

    """
    return (1 - (0.5 * np.exp(-k1*tt) + 0.5 * np.exp(-k*tt))) - A