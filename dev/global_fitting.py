import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import pyhdx
from pyhdx.plot import plot_residue_map, make_coverage_figure
from pyhdx.fitting import KineticsFitting
import numpy as np
import itertools
from symfit import Parameter, Variable, Fit, Model, exp, CallableModel
from symfit.core.minimizers import *
from scipy.optimize import fsolve
from collections import namedtuple
from tqdm.auto import tqdm
from operator import add
from functools import reduce

import pickle

filename = r"C:\Users\jhs\Programming\pyhdx\tests\test_data\ds3.csv"
drop_first = 1  # Number of N terminal residues to ignore

control_100 = ("Full Deuteration control", 0.167)
series_name = 'SecA-monomer'
chisq_thd = 20


pf = pyhdx.PeptideCSVFile(filename, drop_first=drop_first)

# b1 = pf.data['start'] > 97
# b2 = pf.data['start'] < 135
# b = np.logical_and(b1, b2)
# pf.data = pf.data[b]

states = pf.groupby_state_control(control_100)
series = states[series_name]

kf = KineticsFitting(series)
if __name__ == '__main__':
    result = kf.do_fitting(20)

    print('what')

    plt.figure()
    plt.plot(result.rate)
    plt.yscale('log')
    plt.show()


    # with open('temp_seca_full.pick', 'wb') as f:
    #     pickle.dump(result, f)
