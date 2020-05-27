import numpy as np
from pyhdx.simulate import generate_data, trace, to_rates, gen_coverage
from pyhdx import PeptideMasterTable, KineticsFitting
from pyhdx.support import np_from_txt, fmt_export
import matplotlib.pyplot as plt
import pickle

# possible: 7, 15, 24, 27, 29, 31, 38
seed = 10
np.random.seed(seed)


sequence = 'MKIKTGARILALSALTTMMFSASALAKIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLL'


peptides = [
    (5, 15),
    (6, 20),
    (10, 18),
    (12, 30),
]

tr = trace(len(sequence), (15, 10), (50, 100))
rates = to_rates(tr, 10**1.43, 10**-1.2)


timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
cov, data = generate_data(peptides, sequence, timepoints, rates)

np.save('rates.npy', rates)
np.save('full_data.npy', data)

pcf = PeptideMasterTable(data, drop_first=0)
states = pcf.groupby_state()
series = states['state1']


fmt, hdr = fmt_export(data, delimiter=',', width=0)
np.savetxt('gen_data.csv', data, fmt=fmt, header=hdr)

idx = 2

times = [s.exposure for s in series]
upt = [s.data[idx]['scores'] for s in series]

plt.figure()
plt.plot(times, upt)
plt.show()

#
# fitting = KineticsFitting(series)
#
# result1 = fitting.weighted_avg_fit()
# out_arr1 = result1.get_output(['rate', 'tau1', 'tau2', 'r'])
# # np.save(f'fit1_seed_{seed}.npy', out_arr1)
# # with open('fitresult1.pick', 'wb') as f:
# #     pickle.dump(result1, f)
#
# result2 = fitting.lsq_fit_blocks(out_arr1)
# out_arr2 = result2.get_output(['rate', 'tau1', 'tau2', 'r'])
# # np.save(f'fit2_seed_{seed}.npy', out_arr2)
# # with open('fitresult2.pick', 'wb') as f:
# #     pickle.dump(result2, f)
#
# prot_rates = rates[:cov.prot_len]
#
# plt.figure()
# plt.plot(out_arr1['r_number'], prot_rates, label='gt')
# plt.plot(out_arr1['r_number'], out_arr1['rate'], label='fit1')
# plt.plot(out_arr2['r_number'], out_arr2['rate'], label='fit2')
# plt.yscale('log')
# plt.legend()
# plt.show()
