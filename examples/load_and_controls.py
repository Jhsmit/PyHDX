import pyhdx
from pyhdx.fitting import KineticsFitting
from pyhdx.support import np_from_txt
import numpy as np
import pickle
import time
import matplotlib.pyplot as plt
import tensorflow as tf

#tf.Session(config=tf.ConfigProto(allow_growth=True))

filename = r"..\tests\test_data\ds1.csv"
drop_first = 1  # Number of N terminal residues to ignore

control_100 = ('PpiANative', 30.000002)

series_name = 'PpiANative'

data = pyhdx.read_dynamx(filename)
pf = pyhdx.PeptideMasterTable(data, drop_first=drop_first, ignore_prolines=True)
pf.set_control(control_100)
states = pf.groupby_state()
series = states[series_name]
series.make_uniform()
print(series.cov.X.shape)




kf = KineticsFitting(series)
print(series[0].sequence)
#ignore False: 163
#ignore True: 153

#1, True, : (78, 153)

# fit_result = kf.weighted_avg_fit()
#
# t0 = time.time()
# fr1 = kf.weighted_avg_fit()
# t1 = time.time()
#
# print(t1 - t0)
# original: 64s
# TF no seed: 80s (popsize default 40)
# TF seed: 80s (popsize 15)

#arr1 = fr1.get_output(['rate', 'tau', 'k1', 'k2', 'r'])

#np.save('tempresult.npy', arr1)
original = np.load('tempresult.npy')



# with open('fitresult1.pick', 'wb') as f:
#     pickle.dump(fr1, f)

t0 = time.time()
fr2 = kf.blocks_fit(original)
t1 = time.time()
print(fr2)
#original: 200s
#TF no callable model: 180, bounds problems
#original no callable model: > 15 minutes DNF

print(t1 - t0)

arr2 = fr2.get_output(['rate', 'tau', 'k1', 'k2', 'r'])

np.save('fit2.npy', arr2)

fig, ax = plt.subplots()
ax.set_yscale('log')
ax.scatter(arr2['r_number'], arr2['rate'])
ax.scatter(original['r_number'], original['rate'])

plt.show()


# with open('fitresult2.pick', 'wb') as f:
#     pickle.dump(fr2, f)
