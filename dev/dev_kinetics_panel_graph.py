import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import pyhdx
from pyhdx.panel import HDXBase, HDXKinetics
import param
from collections import namedtuple
from threading import Thread
import panel as pn
import numpy as np

pn.extension()



filename = r"C:\Users\jhs\Programming\pyhdx\tests\test_data\ds1.csv"
drop_first = 1  # Number of N terminal residues to ignore
control_state = "PpiA-FD"
control_exposure = .167
series_name = 'PpiANative'

chisq_thd = 20
pf = pyhdx.PeptideCSVFile(filename, drop_first=drop_first)
p_dict = pf.return_by_name(control_state, control_exposure)
p_dict = {k: v for k, v in p_dict.items() if 'Native' in k}

s = p_dict[next(iter(p_dict))]

print(len(s.scores_average))

r_number = np.arange(s.start, s.stop + 1)
print(len(r_number))

hdxk = HDXKinetics(pm_dict=p_dict)


print(hdxk.times)

#
#
hdxk.panel().show(threaded=True)
