import pyhdx
from pyhdx.fitting import KineticsFitting
from pyhdx.support import np_from_txt
import numpy as np

filename = r"..\tests\test_data\ds1.csv"
drop_first = 1  # Number of N terminal residues to ignore

control_100 = ('PpiANative', 30.000002)

series_name = 'PpiANative'

data = pyhdx.read_dynamx(filename)
pf = pyhdx.PeptideCSVFile(data, drop_first=drop_first)
pf.set_control(control_100)
states = pf.groupby_state()
series = states[series_name]


kf = KineticsFitting(series)

fit_result = kf.weighted_avg_fit()
fr1 = kf.weighted_avg_fit()
arr1 = fr1.get_output(['rate', 'tau', 'tau1', 'tau2', 'r'])

fr2 = kf.lsq_fit_blocks(arr1)
arr2 = fr2.get_output(['rate', 'tau', 'tau1', 'tau2', 'r'])
