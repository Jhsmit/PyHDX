import tensorflow as tf
from tensorflow.keras.layers import Layer
from tensorflow.keras.constraints import Constraint
from tensorflow.keras.regularizers import Regularizer
from tensorflow.python.keras import backend as K
from tensorflow.python.ops import math_ops


class L1L2Differential(Regularizer):
    def __init__(self, l1=0., l2=0.):  # pylint: disable=redefined-outer-name
        self.l1 = K.cast_to_floatx(l1)
        self.l2 = K.cast_to_floatx(l2)

    def __call__(self, x):
        if not self.l1 and not self.l2:
            return K.constant(0.)
        regularization = 0.
        diff = x[:-1] - x[1:]
        if self.l1:
            regularization += self.d*math_ops.reduce_sum(math_ops.abs(diff))
        if self.l2:
            regularization += self.d * math_ops.reduce_sum(math_ops.square(diff))

        return regularization

    def get_config(self):
        return {'l1': float(self.l1), 'l2': float(self.l2)}


class Between(Constraint):
    def __init__(self, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value

    def __call__(self, w):
        return tf.clip_by_value(w, self.min_value, self.max_value)

    def get_config(self):
        return {'min_value': self.min_value,
                'max_value': self.max_value}


class TFParameter(object):
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

    def __init__(self, timepoints, kint):
        self.timepoints = tf.dtypes.cast(tf.expand_dims(timepoints, 0), tf.float32)
        self.kint = tf.dtypes.cast(tf.expand_dims(kint, -1), tf.float32)

    def __call__(self, X, **parameters):
        pfact = 10**parameters['log_P']
        uptake = 1 - tf.exp(-tf.matmul((self.kint/pfact), self.timepoints))
        return 100*tf.matmul(X, uptake)

    def compute_output_shape(self, input_shape):
        return input_shape[0], len(self.timepoints)


class AssociationRateFunc(object):

    def __init__(self, timepoints):
        self.timepoints = tf.dtypes.cast(tf.expand_dims(timepoints, 0), tf.float32)

    def __call__(self, X, **parameters):
        k = 10**parameters['k']
        uptake = 1 - tf.exp(-tf.matmul(k, self.timepoints))
        return 100*tf.matmul(X, uptake)

    def compute_output_shape(self, input_shape):
        return input_shape[0], len(self.timepoints)



class TFFitResult(object):
    """


    Parameters
    ----------
    r_number list or r numbers these results cover
    intervals (inclusive, exclusive) intervals which map results, models to r numbers
    results list of results returned from model.fit
    models list of tensorflow models
    """
    def __init__(self, r_number, intervals, results, models):
        assert len(results) == len(models)
        #        assert len(models) == len(block_length)
        self.r_number = r_number
        self.intervals = intervals  # inclusive, excluive
        self.results = results
        self.models = models


