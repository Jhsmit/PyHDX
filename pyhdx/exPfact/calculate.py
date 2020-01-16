import numpy as np

from scipy import optimize

from calc_dpred import calculate_dpred


def calculate_rms(dpred, dexp, nj, weights=None):
    """
    Calculates the normalised rms between dpred and dexp
    :param dpred: numpy array containing the dpred values.
    :param dexp: numpy array containing the dexp values.
    :param nj: number of residues in the protein.
    :return: rms (float)
    """
    if weights is not None:
        return 1 / nj * (np.sum([weights[i]*(dpred[i] - dexp[i]) ** 2 for i in range(len(dexp))]))
    else:
        return 1 / nj * (np.sum([(dpred[i] - dexp[i]) ** 2 for i in range(len(dexp))]))


def cost_function(params, *args):
    """
    Cost function for pfactor fitting.
    :param params: list of estimated pfactors
    :param args: arguments required for calculating the cost
    :return: cost score (float)
    """
    dexp, tk, assignments, k, kint, weights = args
    if weights is not None:
        score = calculate_rms(calculate_dpred(np.array(params), tk, kint, assignments), dexp, len(assignments), weights)
    else:
        score = calculate_rms(calculate_dpred(np.array(params), tk, kint, assignments), dexp, len(assignments))
    score += harmonic_score(params, k)
    return float(score)


def do_random_search(kint, search_steps, pfactor_filter, dexp, time_points, assignments, harmonic_term, prolines, weights):
    """

    :param kint: array of kint values.
    :param search_steps: integer of steps to take.
    :param assignment_set: Python Set of assignments.
    :param dexp: array of dexp values.
    :param time_points: array of time points.
    :param assignments: array of assignment arrays.
    :param harmonic_term: term to use for harmonic scoring.
    :return: dictionary containing all scores mapped to pfactor arrays.
    """

    score_array = {}
    for i in range(search_steps):

        init_array = [
            np.random.uniform(0.01, 20.00) if ii != 0 and ii+1 not in prolines and ii+1 in pfactor_filter else -1
            for ii in range(max(pfactor_filter))]

        score = cost_function(init_array, dexp, time_points,
                              assignments, harmonic_term, kint, weights)

        #my run:
        # init_array: list
        # dexp:tuple with arrays (len(2)
        # time_points: array

        score_array[score] = init_array

    return score_array


def fit_pfact(init_array, dexp, time_points, assignments, harmonic_term, kint, bounds, tol, weights):
    """
    :param init_array: initial guess array of pfactors for minimisation.
    :param dexp: array of dexp values.
    :param time_points: array of time points.
    :param kint: array of kint values.
    :param assignments: array of assignment arrays.
    :param harmonic_term: term to for harmonic cost scoring.
    :param bounds: array of paired bounds for minimisation.
    :param tol: tolerance for minimisation convergence
    :return: scipy assignment object containing optimum
    """
    pfit = optimize.minimize(cost_function, init_array,
                             args=(dexp,
                                   time_points,
                                   assignments,
                                   harmonic_term,
                                   kint,
                                   weights),
                             method='L-BFGS-B',
                             bounds=bounds,
                             tol=tol,
                             options={'disp': True,
                                      'maxfun': 150000000,
                                      'maxiter': 1500000
                                      }
                             )
    return pfit


def harmonic_score(params, k):
    """
    Calculates harmonic score for pfactors.
    :param params: array of pfactors.
    :param k: harmonic term.
    :return: score (float)
    """
    scores = []
    for ii in range(len(params) - 1):
        if params[ii] >= 0 and params[ii + 1] >= 0:
            scores.append((k / 2) * (params[ii] - params[ii + 1]) ** 2)
    return sum(scores)


def predict_dexp(pfact, time_points, assignments):
    """
    Calculates predicted dexp from pfactors, time_points, assignments and kint values.
    :param pfact: array of pfactors.
    :param time_points: array of time points.
    :param kint: array of kint values
    :param assignments: array of assignment arrays.
    :return: numpy array of dexp values.
    """
    dexp = calculate_dpred(pfact, time_points, assignments)
    return dexp
