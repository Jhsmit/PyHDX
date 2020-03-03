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
b = data['end'] < 150
data = data[b]
pf = pyhdx.PeptideCSVFile(data, drop_first=drop_first)
pf.set_control(control_100)
states = pf.groupby_state()
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
print(len(sec_res))
kf_sec = KineticsFitting(sec1)
print(kf_sec.scores_stack.shape)
print(kf_sec.scores_stack)

print('howdoe')
#
# print(sec1.cov.block_length)
# print(sec_res)

#
# class LSQKinetics(KineticsModel): #TODO find a better name (lstsq)
#     #todo block length is redundant
#     def __init__(self, initial_result, kf_section):
#         """
#
#         Parameters
#         ----------
#         initial_result kineticsresult object from initial fitting
#         kf_section kineticsfitting object for the section
#         """
#         super(LSQKinetics, self).__init__()
#         # print(block_lengths)
#         # total = np.sum(block_lengths)
#         t_var = self.make_variable('t')
#
#         #Assemble terms spanning accross blocks of residues
#         #blocks with only 1 residue have only 1 time component
#         terms = []
#         for i, (r, m, bl) in enumerate(initial_result):
#             if bl == 1:
#                 t1v = r.params[m.names['tau1']]
#                 t2v = r.params[m.names['tau2']]
#                 rv = r.params[m.names['r']]
#                 value = rv * t1v + (1 - rv) * t2v
#                 tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1 / 40, value=value)
#
#                 term = (1 - exp(-t_var / tau1))
#             else:
#                 t1v = r.params[m.names['tau1']]
#                 t2v = r.params[m.names['tau2']]
#                 rv = r.params[m.names['r']]
#                 tau1 = self.make_parameter('tau1_{}'.format(i), max=30, min=1 / 40, value=t1v)
#                 tau2 = self.make_parameter('tau2_{}'.format(i), max=30, min=1 / 40, value=t2v)
#                 r = self.make_parameter('r_{}'.format(i), max=1, min=0, value=rv)
#
#                 term = (1 - (r*exp(-t_var / tau1) + (1-r)*exp(-t_var/tau2)))
#             terms.append(term)
#
#         #Iterate over rows (peptides) and add terms together which make one peptide
#         model_dict = {}
#         d_vars = []
#         for i, x_row in enumerate(kf_section.k_series.cov.X_red_norm):
#             d_var = self.make_variable('d_{}'.format(i))
#             d_vars.append(d_var)
#             rhs = reduce(add, [100*fraction*term for fraction, term in zip(x_row, terms)])
#             model_dict[d_var] = rhs
#
#         self.d_vars = d_vars
#         self.t_var = t_var
#         self.model = CallableModel(model_dict)
#

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

                component = 100 * (block / total) * (r*exp(-t_var / tau1) + (1-r)*exp(-t_var/tau2) )

            print(block, component)
            d_var = self.make_variable('d_{}'.format(i))
            model_dict[d_var] = component

        self.model = CallableModel(model_dict)


make_coverage_figure(sec1.cov, 5, 75)
#plt.show()
plt.savefig('testcoverage.png')

m = LSQKinetics(sec_res, kf_sec)

plt.figure()
plt.plot(sec_res.rate)
plt.show()
print(m.sf_model)

sp = kf_sec.scores_peptides
print(sp.shape)

data_dict = {d_var.name: scores for d_var, scores in zip(m.d_vars, kf_sec.scores_peptides.T)}
print(data_dict)
data_dict[m.t_var.name] = sec1.times
print('halla')

fit = Fit(m.sf_model, **data_dict)
initial_guesses = {par.name: par.value for par in m.sf_model.params}
initial_guesses[m.t_var.name] = sec1.times

output = m.sf_model(**initial_guesses).output_dict

initial_chisq = reduce(add, [np.sum((v - data_dict[k.name])**2) for k, v in output.items()])
print(initial_chisq)

fig, ax = plt.subplots()
for k, v in output.items():
    ax.plot(sec1.times, v, label=k)
    ax.scatter(sec1.times, data_dict[k.name])
plt.xscale('log')
plt.legend()
plt.xlabel('Time (min)')
plt.ylabel('Deuteration (%)')
plt.show()

res = fit.execute()
print(res)

output = m.sf_model(**{m.t_var.name: sec1.times}, **res.params).output_dict

final_chisq = reduce(add, [np.sum((v - data_dict[k.name])**2) for k, v in output.items()])
print(final_chisq)


fig, ax = plt.subplots()
for k, v in output.items():
    ax.plot(sec1.times, v, label=k)
    ax.scatter(sec1.times, data_dict[k.name])
plt.xscale('log')
plt.legend()

plt.xlabel('Time (min)')
plt.ylabel('Deuteration (%)')
plt.show()

for k, v in res.params.items():
    print(k, v)

rate = m.get_rate(**res.params)

r_number = sec1.cov.r_number
i0, i1 = np.searchsorted(r_number, [s, e])
print(r_number[i0:i1+1])

r_number[i0:i1+1] = rate

plt.figure()
plt.plot(rate)
plt.show()



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
