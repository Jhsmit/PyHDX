import numpy as np
from pyhdx.simulate import generate_data, trace, to_rates, gen_coverage
from pyhdx import PeptideMasterTable, KineticsFitting, read_dynamx
import pprint as pp
from pyhdx.support import np_from_txt, fmt_export
import matplotlib.pyplot as plt
import pickle


# possible: 7, 15, 24, 27, 29, 31, 38
seed = 10
np.random.seed(seed)


sequence = 'MKIKPTPPRILALSAPLTTMMFSASALAPKIEEGKLVIPWINGDKGYNPGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLL'

peptides = [
    (5, 15),   # length 11 (inclusive, inclusive)
    (6, 20),  #15
    (10, 18), #9
    (12, 30),
    (35, 42),
    (37, 45),
    (40, 44)
]
#[11 15  9 19  8  9  5]
#
# [( 5, 15,   0.167, 'state1', 'TPPRILALSAP',   4.1213192)
#  ( 6, 20,   0.167, 'state1', 'PPRILALSAPLTTMM',   4.0046056)
#  (10, 18,   0.167, 'state1', 'LALSAPLTT',   4.1213192)
#  (12, 30,   0.167, 'state1', 'LSAPLTTMMFSASALAPKI',   3.1077533)
#  (35, 42,   0.167, 'state1', 'LVIPWING',  63.38898  )
#  (37, 45,   0.167, 'state1', 'IPWINGDKG',  89.15229  )
#  (40, 44,   0.167, 'state1', 'INGDK', 100.       )

tr = trace(len(sequence), (15, 10), (50, 100))
rates = to_rates(tr, 10**1.43, 10**-1.2)


timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
cov, data = generate_data(peptides, sequence, timepoints, rates)

print(data.sequence)


raise

print(data.dtype)

np.save('rates.npy', rates)
np.save('full_data.npy', data)

fmt, hdr = fmt_export(data, delimiter=',', width=0)
print(hdr)
np.savetxt('simulated_data.csv', data, fmt=fmt, header=hdr, delimiter=',')
np.savetxt('simulated_rates.txt', rates)
