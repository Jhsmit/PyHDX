import tensorflow as tf
from tensorflow.keras.layers import Layer, Input
from tensorflow.keras.constraints import Constraint
from tensorflow.keras.regularizers import Regularizer
from tensorflow.keras.initializers import Initializer
from tensorflow.python.keras import backend as K
from tensorflow.python.ops import math_ops
from tensorflow.keras.optimizers import Adagrad
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.losses import Loss
from pyhdx.models import Protein
import numpy as np
import copy


class NaNMeanSquaredError(Loss):
    """MSE which ignores nan entries

    Parameters
    ----------
    y_true:

    """
    def call(self, y_true, y_pred):
        y_pred = tf.convert_to_tensor(y_pred)
        y_true = math_ops.cast(y_true, y_pred.dtype)
        diff = math_ops.square(y_pred - y_true)  # has nans

        nanmean = tf.reduce_mean(tf.boolean_mask(diff, ~tf.math.is_nan(diff)))
        return nanmean


class L1L2Differential(Regularizer):
    """
    A regularized that applies and L1 or L2 regularization penalty to the differential of a parameter vector.

    Parameters
    ----------
    l1: :obj:`float`
        L1 regularization factor
    l2: :obj:`float`
        L2 regularization factor

    """
    def __init__(self, l1=0., l2=0.):  # pylint: disable=redefined-outer-name
        self.l1 = K.cast_to_floatx(l1)
        self.l2 = K.cast_to_floatx(l2)

    def __call__(self, x):
        if not self.l1 and not self.l2:
            return K.constant(0.)
        regularization = 0.
        diff = x[:-1] - x[1:]
        if self.l1:
            regularization += self.l1*tf.math.reduce_mean(math_ops.abs(diff))
        if self.l2:
            regularization += self.l2*tf.math.reduce_mean(math_ops.square(diff))  # should we take the sqrt of the sum?

        return regularization

    def get_config(self):
        return {'l1': float(self.l1), 'l2': float(self.l2)}


class Between(Constraint):
    """
    Interval parameter constraint.

    Constrains the values of parameters to the interval [min_value, max_value].

    Parameters
    ----------
    min_value: :obj:`float`
        Lower bound for the allowed interval (optional `None`).
    max_value: :obj:`float`
        Upper bound for the allowed interval (optional `None`).
    """
    def __init__(self, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value

    def __call__(self, w):
        return tf.clip_by_value(w, self.min_value, self.max_value)

    def get_config(self):
        return {'min_value': self.min_value,
                'max_value': self.max_value}


class TFParameter(object):
    """
    Parameter objects used in `CurveFit` TensorFlow Layer.
    Parameters are 'weights' in the context of Neural Networks.

    Parameters
    ----------
    name: :obj:`str`
        Name of the parameter
    shape: :obj:`tuple`
        Parameter shape
    initializer: :class:`~tensorflow.python.keras.initializers.Initializer`
        Subclass of Keras Initializer to initialize parameter elements.
    regularizer :class:`~tensorflow.python.keras.regularizers.Regularizer`
        Subclass of Keras Regularizer applied to parameter elements.
    constraint :class:`~tensorflow.python.keras.constraints.Constraint`
        Subclass of keras Constraint applied to parameter elements.

    """
    def __init__(self, name, shape, initializer=None, regularizer=None, constraint=None):
        self.name = name
        self.shape = shape
        self.initializer = initializer
        self.regularizer = regularizer
        self.constraint = constraint


class CurveFit(Layer):
    def __init__(self, params, function, **kwargs):
        super(CurveFit, self).__init__(**kwargs)
        self.params = params
        self.function = function

        self.parameters = {}

    def build(self, input_shape):
        for param in self.params:
            wts = self.add_weight(name=param.name,
                                  shape=param.shape,
                                  regularizer=param.regularizer,
                                  initializer=param.initializer,
                                  constraint=param.constraint,
                                  trainable=True)
            self.parameters[param.name] = wts

        super().build(input_shape)

    def call(self, inputs, **kwargs):
        return self.function(inputs, **self.parameters)

    def compute_output_shape(self, input_shape):
        return self.function.compute_output_shape(input_shape)


class LossHistory(tf.keras.callbacks.Callback):

    def __init__(self, verbose=False):
        super(LossHistory, self).__init__()
        self.verbose = verbose
        self.current_loss = 1.e20
        self.weights = None

    def on_epoch_end(self, epoch, logs=None):
        loss = logs.get('loss')
        if logs.get('loss') < self.current_loss:
            self.current_loss = loss
            self.weights = list([layer.get_weights() for layer in self.model.layers])

        if self.verbose and epoch % 1000 == 0:
            print(epoch, loss)


class AssociationPFactFunc(object):
    parameter_name = 'log_P'

    def __init__(self, timepoints):
        self.timepoints = tf.dtypes.cast(tf.expand_dims(timepoints, 0), tf.float32)

    def __call__(self, inputs, **parameters):
        """

        Parameters
        ----------
        inputs: :obj:`list`
        list of X, k_int, shapes (N_peptides, N_residues), (N_peptides, 1)
        parameters: dict with fit parameters

        Returns
        -------

        """
        pfact = 10**parameters[self.parameter_name]
        uptake = 1 - tf.exp(-tf.matmul((inputs[1]/(1 + pfact)), self.timepoints))
        return tf.matmul(inputs[0], uptake)

    # def compute_output_shape(self, input_shape):
    #     return input_shape[0], len(self.timepoints)

    def call_numpy(self, inputs, **parameters):
        pfact = 10**parameters[self.parameter_name]
        uptake = 1 - np.exp(-np.matmul((inputs[1]/(1 + pfact)), np.array(self.timepoints)))
        return np.matmul(inputs[0], uptake)

    @staticmethod
    def output(weights):
        #todo perhaps input weights as dict (and store as such in HistoryCallback)
        return {'pfact', 10**weights}


class AssociationRateFunc(object):
    parameter_name = 'log_k'

    """
    Function passed to CurveFit layer to calculate forward propagation through the network



    """
    def __init__(self, timepoints):

        self.timepoints = tf.dtypes.cast(tf.expand_dims(timepoints, 0), tf.float32)

    def __call__(self, inputs, **parameters):
        k = 10**parameters[self.parameter_name]
        uptake = 1 - tf.exp(-tf.matmul(k, self.timepoints))
        return tf.matmul(inputs[0], uptake)

    # def compute_output_shape(self, input_shape):
    #     return input_shape[0], len(self.timepoints)


class TFFitResult(object):
    """

    Parameters
    ----------
    r_number list or r numbers these results cover
    intervals (inclusive, exclusive) intervals which map results, models to r numbers  (can be obtained from series)
    funcs: assumed to be tghe same
    assumed to be the same for all intervals
    weights: list of weights (parameters) at lowest loss
    """
    def __init__(self, series, intervals, funcs, weights, inputs, loss=None):
        #todo remove intervals
        self.r_number = series.cov.r_number
        self.series = series
        self.intervals = intervals  # inclusive, excluive
        self.funcs = funcs
        #self.results = results
        self.weights = weights
        self.inputs = inputs

        self.loss = loss

    @property
    def output(self):
        #todo merge with KineticsFitresults get_param function

        name = self.funcs[0].parameter_name
        output = np.full_like(self.r_number, fill_value=np.nan, dtype=float)

        for (s, e), wts in zip(self.intervals, self.weights):  #there is only one interval in the current implementation
            i0, i1 = np.searchsorted(self.r_number, [s, e])
            output[i0:i1] = wts

        array = np.empty_like(self.r_number, dtype=[('r_number', int), (f'{name}_full', float), (name, float)])
        array['r_number'] = self.r_number
        array[f'{name}_full'] = output

        bools = ~self.series.cov['exchanges']
        array[name] = output.copy()
        array[name][bools] = np.nan  # set no coverage or prolines resiudes to nan

        return Protein(array, index='r_number')

    def __call__(self, timepoints):
        """output: N x M array (peptides, timepoints)"""
        d_list = []
        for func, wts, ip in zip(self.funcs, self.weights, self.inputs):
            f_copy = copy.copy(func)
            f_copy.timepoints = tf.dtypes.cast(tf.expand_dims(timepoints, 0), tf.float32)

            parameters = {func.parameter_name: tf.dtypes.cast(tf.expand_dims(wts, -1), tf.float32)}
            X = tf.dtypes.cast(np.squeeze(ip[0], axis=0), tf.float32)
            k_int = tf.dtypes.cast(np.squeeze(ip[1], axis=0), tf.float32)
            d_out = f_copy([X, k_int], **parameters)

            d = np.array(d_out)
            d_list.append(d)

        full_output = np.concatenate(d_list)

        return full_output


