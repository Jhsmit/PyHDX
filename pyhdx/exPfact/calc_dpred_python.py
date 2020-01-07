from math import exp
import numpy as np


def get_rate_res(assignments, pfactors):


    kint_list = np.zeros((len(assignments), 50))
    for assignment_index in range(len(assignments)):

        start_index = int(assignments[assignment_index][1])
        end_index = int(assignments[assignment_index][2])
        N = end_index - start_index
        pfactor_subset = pfactors[start_index:end_index + 1]
        kint_list[assignment_index][0] = assignments[assignment_index][0]
        kint_list[assignment_index][1] = assignments[assignment_index][1]
        kint_list[assignment_index][2] = assignments[assignment_index][2]
        for i in range(3, N):
            if assignments[assignment_index][i] == -1:
                kint_list[assignment_index][i] += -1
            else:
                kint_list[assignment_index][i] += 60 * assignments[assignment_index][i + 3] / exp(pfactor_subset[i])

    return np.array(kint_list)


def res2frag(kint_list, time_points):

    nfrag = int(kint_list.shape[0])
    ntime = int(time_points.shape[0])
    rate_frag = np.zeros((nfrag, ntime))
    namide = np.zeros(nfrag)

    for i in range(nfrag):
        value_range = int(kint_list[i][2] - kint_list[i][1])
        for j in range(value_range):
            if (kint_list[i][j + 3] >= 0):
                namide[i] = namide[i] + 1

    for i in range(nfrag):
        for k in range(ntime):
            value_range = int(kint_list[i][2] - kint_list[i][1])
            for j in range(3, value_range):
                if (kint_list[i][j + 3] >= 0):
                    rate_frag[i][k] = rate_frag[i][k] + exp(-kint_list[i][j + 3] * time_points[k])
            rate_frag[i][k] = (namide[i] - rate_frag[i][k]) / namide[i]

    return rate_frag


def calculate_dpred(pfactors,
                    time_points,
                    assignments):

    kint_list = get_rate_res(assignments, pfactors)

    frag_rates = res2frag(kint_list, time_points)

    return frag_rates