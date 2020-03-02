import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import pyhdx
from pyhdx.plot import plot_residue_map, make_coverage_figure
from pyhdx.fitting import KineticsFitting, KineticsModel
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

control_100 = ("Full Deuteration control", 0.167)
series_name = 'SecA-monomer'
chisq_thd = 20

data = pyhdx.read_dynamx(filename)
b = data['end'] < 150
data = data[b]
pf = pyhdx.PeptideCSVFile(data, drop_first=drop_first)

states = pf.groupby_state_control(control_100)
series = states[series_name]
series.make_uniform()
split = series.split()
sec1 = split['31_50']


total_cs = np.append(series.cov.start, series.cov.start + np.cumsum(series.cov.block_length))
i0, i1 = np.searchsorted(total_cs, [31, 50])

pm0 = sec1[0]



#print(series[0].coverage)
print('cov len', len(series[0].has_coverage))
kf = KineticsFitting(series)



# if __name__ == '__main__':
#     kf = KineticsFitting(series)
#     fr = kf.do_fitting()
#     print(len(fr.results))
#
#
#     with open('temp_fit.pick', 'wb') as f:
#         pickle.dump(fr, f)
#
with open('temp_fit.pick', 'rb') as f:
    fit_res = pickle.load(f)

print(total_cs)
#
print(split.keys())
sec1 = split['31_50']
s, e = sec1.cov.start, sec1.cov.end
i0, i1 = np.searchsorted(total_cs, [s, e + 1])
#sec_res = fit_res.results[i0:i1]
sec_res = fit_res[i0:i1]
print(sec_res)


# print(sec_res)
#
# print(sec1.cov.block_length)
# print(sec_res)
# #
class OverallKinetics(KineticsModel): #TODO find a better name (lstsq)
    #todo block length is redundant
    def __init__(self, initial_result, block_lengths):
        super(OverallKinetics, self).__init__()
        assert len(initial_result) == len(block_lengths)
        print(block_lengths)
        total = np.sum(block_lengths)
        model_dict = {}

        t_var = self.make_variable('t')
        for i, (res, block) in enumerate(zip(initial_result, block_lengths)):
            r, m, b = res
            if block == 1:
                #Only one time component for residue blocks of 1
                t1v = r.params[m.names['tau1']]
                t2v = r.params[m.names['tau2']]
                rv = r.params[m.names['r']]
                value = rv*t1v + (1-rv)*t2v

                tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1/40, value=value)

                component = 100 * (block / total) * exp(-t_var / tau1)
            else:
                t1v = r.params[m.names['tau1']]
                t2v = r.params[m.names['tau2']]
                rv = r.params[m.names['r']]
                tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1 / 40, value=t1v)
                tau2 = self.make_parameter('tau2_{}'.format(i), max=30, min=1 / 40, value=t2v)
                r = self.make_parameter('r_{}'.format(i), max=1, min=0, value=rv)

                component = 100 * (block / total) * ( r*exp(-t_var / tau1) + (1-r)*exp(-t_var/tau2) )

            print(block, component)
            d_var = self.make_variable('d_{}'.format(i))
            model_dict[d_var] = component

        self.model = CallableModel(model_dict)


m = OverallKinetics(sec_res, sec1.cov.block_length)
print(m.model)

# #assert np.all(sec1.cov.has_coverage == series.cov.has_coverage)
# # for section in split.values():
# #     s, e = section.cov.start, section.cov.end
# #     i0, i1 = np.searchsorted(total_cs, [s, e + 1])
# #     sec_res = initial_results.results[i0:i1]
# # #
# #     print(sec_res)
# # #
# #     for res in sec_res:
# #         try:
# #             assert isinstance(res, FitResults)
# #         except AssertionError:
# #             print(s, e)
# #             start = series.cov.start
# #             print(sec_res)
#
#     # with open('temp_seca_full.pick', 'wb') as f:
#     #     pickle.dump(result, f)
