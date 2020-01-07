from read import read_assignments, read_kint, read_pfact, read_time_points

from calculate import predict_dexp

import sys

import numpy as np

import re

import os


def get_dcalc(assignments_file, kint_file, pfact_file, time_points_file, fragment_number):

    assignments = read_assignments(assignments_file)
    kint = read_kint(kint_file, -1)
    pfact = read_pfact(pfact_file)
    time_points = read_time_points(time_points_file)

    residue_range = np.array([assignments[int(fragment_number)-1]])[0]

    dcalc_all = np.zeros((len(time_points), residue_range[2]-residue_range[1]+1))

    for ii, time in enumerate(time_points):
        k = kint[residue_range[1]:residue_range[2]]
        p = pfact[residue_range[1]:residue_range[2]]
        dcalc = np.insert(1.0 - np.exp(-k * 60 * time / p), 0, time)
        dcalc_all[ii] = dcalc

    return dcalc_all

if __name__ == "__main__":

    assignments_file, kint_file, pfact_file, time_points_file, fragment_number = sys.argv[1:]
    dcalc = get_dcalc(assignments_file, kint_file, pfact_file, time_points_file, fragment_number)
    n = dcalc.shape[1]-1

    for ii, val in enumerate(dcalc):

        os.system("sed -e 's/KKK/{}/g' ../R/binomial.tmpl > r.inp".format(','.join([str(x) for x in val[1:].tolist()])))
        os.system("R --no-save < r.inp > r.out")

        # iso_file = open("{}.iso".format(fragment_number), 'w')
        iso_lines = []
        with open("{}.prospector.ucsf.edu".format(fragment_number)) as f:

            for line in f.readlines():
                if not line.startswith('#') and float(line.split()[2]) > 0:
                    ll = line.split()
                    iso_line = (int(ll[0]), float(ll[2])/100, float(ll[1]))
                    iso_lines.append(iso_line)

        with open('r.out', 'r') as f:
            for line in f.readlines():
                if line.startswith('Pr'):
                    l = line.split()
                    iso_lines.append((int(l[1]), float(l[2])))


        two_values = [x for x in iso_lines if len(x) == 2]
        three_values = [x for x in iso_lines if len(x) == 3]


        probs = [("{:.2f}".format(val1[0] * 1.00627 + val2[2]), val1[1] * val2[1]) for val1 in two_values for val2 in three_values]

        probs.sort(key=lambda x: x[0])

        fout = open("{}.ist".format(ii+1), 'w')

        for x in sorted(set([x[0] for x in probs])):
            y_array = []
            for y in probs:
                if y[0] == x:
                    y_array.append(y[1])
            fout.write(' '.join([str(x), str(sum(y_array)), '\n']))
        fout.close()