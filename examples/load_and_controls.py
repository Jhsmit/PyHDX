import pyhdx
from pyhdx.fitting import KineticsFitting

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
