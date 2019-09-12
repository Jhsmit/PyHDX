import numpy as np
import scipy.optimize


# from https://github.com/Fluorescence-Tools/ChiSurf/blob/master/mfm/math/optimization/nnls.py
def solve_nnls(A, b, reg=1e-5, x_max=1e12):
    """
    Solve :math:`argmin_x || Ax - b ||_2 for x>=0`. This is a wrapper for a FORTAN non-negative least squares solver (as
    in the SciPy. In addition to that the matrix A is regularized by a Thikonov-regularization. If values bigger
    than x_max are encountered in the solution they are simply replaced by zeros.
    :param A: numpy-array
        Matrix A as shown above.
    :param b: numpy-array
        Right-hand side vector.
    :param reg: float
        Regularization factor
    :param x_max: float
    :return:
    """
    right = np.dot(A, b)
    left = np.dot(A, A.T) + reg * np.diag(np.ones_like(right)) * A.T.shape[0] / A.T.shape[1]
    x = scipy.optimize.nnls(left, right)[0]
    x[x > x_max] = 0
    return x
