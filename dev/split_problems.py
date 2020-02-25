import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import pyhdx
from pyhdx import KineticsSeries
from pyhdx.plot import plot_residue_map, make_coverage_figure
from pyhdx.fitting import KineticsFitting
from pyhdx.support import reduce_inter
from functools import reduce
from operator import add
import numpy as np

import pickle

filename = r"C:\Users\jhs\Programming\pyhdx\tests\test_data\ds3.csv"
drop_first = 1  # Number of N terminal residues to ignore

control_100 = ("Full Deuteration control", 0.167)
series_name = 'SecA-monomer'
chisq_thd = 20

pf = pyhdx.PeptideCSVFile(filename, drop_first=drop_first)

states = pf.groupby_state_control(control_100)
series = states[series_name]

pm = series[0]
small = [(s, e) for s, e in zip(pm.data['start'], pm.data['end'])]

r1 = reduce_inter(small)

print(len(pm.data))

full_ds = np.concatenate([pm.data for pm in series])
print(len(full_ds))

output = {}
for section in r1:
    s, e = section
    b = np.logical_and(full_ds['start'] >= s, full_ds['end'] <= e)
    output['{}_{}'.format(s, e)] = KineticsSeries(full_ds[b])


new_len = reduce(add, [reduce(add, [len(pm.data) for pm in ks]) for ks in output.values()])
print(new_len)


sets = [{(s, e, seq) for s, e, seq in zip(pm.data['start'], pm.data['end'], pm.data['sequence'])} for pm in series]
print(sets)

intersection = set.intersection(*sets)

print(intersection)
print(len(intersection))

dtype = series.full_data['sequence'].dtype

arr = np.array([tup for tup in intersection], dtype=[('start', int), ('end', int), ('sequence', series.full_data['sequence'].dtype)])
print(arr)
print(len(arr))
#
# b = arr == pm.data[['start', 'end', 'sequence']]
# print(b)


#
def fields_view(arr, fields):
    dtype2 = np.dtype({name: arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)


vf = fields_view(series[0].data, ['start', 'end', 'sequence'])
vf1 = fields_view(arr, ['start', 'end', 'sequence'])

isin = np.isin(vf, vf1)
print(len(vf), len(vf1))
isin = np.isin(series[0].data[['start', 'end', 'sequence']], arr)
print(isin)

#print(vf == vf1)
print(pm.data.dtype.fields)
print([pm.data.dtype.fields[name] for name in ['start', 'end', 'sequence']])

print(series.uniform)



series.make_uniform()
print(series.uniform)

for pm in series:
    print(len(pm))


