from pyhdx.fileIO import read_dynamx, txt_to_np, fmt_export
from pyhdx.models import PeptideMasterTable
from numpy.lib.recfunctions import append_fields
import numpy as np

array = txt_to_np('test_data/simulated_data.csv', delimiter=',')
data = read_dynamx('test_data/simulated_data.csv')
print(array.dtype.names)
print(data.dtype.names)

pmt = PeptideMasterTable(array, drop_first=0, ignore_prolines=False, remove_nan=False)

print(pmt.data.dtype.names)

uptake = pmt.data['ex_residues'] * pmt.data['scores'] / 100

for u, m in zip(uptake, pmt.data['ex_residues']):
    print(u, m)

extended = append_fields(pmt.data, ['uptake'], [uptake], usemask=False)
fields = ('start', 'end', 'exposure', 'state', 'sequence', 'ex_residues', 'uptake')

dtype = [(name, extended[name].dtype) for name in fields]
export = np.empty_like(uptake, dtype=dtype)
for name in fields:
    export[name] = extended[name]
fmt, hdr = fmt_export(export, delimiter=',', width=0)
np.savetxt('test.txt', export, fmt=fmt, header=hdr)


new_data = read_dynamx('test.txt')


pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)
states = pmt.groupby_state()
series = states['state1']
series.make_uniform()

