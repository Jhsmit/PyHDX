import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import pyhdx
from pyhdx.plot import plot_residue_map, make_coverage_figure
from pyhdx.fitting import KineticsFitting, KineticsModel, LSQKinetics
import numpy as np
import itertools
from symfit import Parameter, Variable, Fit, Model, exp, CallableModel, FitResults
from symfit.core.minimizers import *
from scipy.optimize import fsolve
from collections import namedtuple
from tqdm.auto import tqdm
from operator import add
from functools import reduce

import pickle

filename = r"C:\Users\jhs\Programming\pyhdx\tests\test_data\ds3.csv"
drop_first = 1  # Number of N terminal residues to ignore

#control_100 = ("Full Deuteration control", 0.167)
control_100 = ("SecA-monomer", 30.000002)
series_name = 'SecA-monomer'
chisq_thd = 20

data = pyhdx.read_dynamx(filename)
# b = data['end'] < 150
# data = data[b]
pf = pyhdx.PeptideCSVFile(data, drop_first=drop_first)
pf.set_control(control_100)
states = pf.groupby_state()
series = states[series_name]
series.make_uniform()
#
# if __name__ == '__main__':
#     kf = KineticsFitting(series)
#     fr = kf.do_fitting()
#     print(len(fr.results))
#
#
#     with open('full_fit_secA_0303.pick', 'wb') as f:
#         pickle.dump(fr, f)


with open('full_fit_secA_0303.pick', 'rb') as f:
    fit_res = pickle.load(f)
#
#
kf_overall = KineticsFitting(series)
r, rate = kf_overall.fine_fitting(fit_res)


plt.figure()
plt.plot(r, rate)
plt.yscale('log')
plt.show()

