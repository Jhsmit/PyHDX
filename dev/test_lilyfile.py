from pyhdx import PeptideCSVFile
import numpy as np

#fpath = os.path.join(directory, 'test_data', 'ds2.csv')
filename = r"C:\Users\jhs\Programming\pyhdx\tests\test_data\ds2.csv"

pf2 = PeptideCSVFile(filename)
filename = "wt_ppiA_folding_4Cmodif_230120.csv"
drop_first = 1  # Number of N terminal residues to ignore
control_state = "FD"
control_exposure = 0.001
series_name = 'folding_4C_10secLabelling'

control_100 = ("FD", 0.001)
control_0 = ("Native folded", 60.000004)

chisq_thd = 20


states = pf2.groupby_state_control(control_100, remove_nan=False)
print(states.keys())
series = states['folding_4C_10secLabelling']

print(series[1].scores)
series[1].scores[:] = 100
print(series[1].scores)

print(series[1].scores_average)
pm = series[1]
avg = pm.X_norm.T.dot(pm.scores)

print(avg)
