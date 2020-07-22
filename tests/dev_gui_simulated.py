from pyhdx.panel.controller import Controller
from pyhdx.fileIO import read_dynamx
from pyhdx import PeptideMasterTable
from pyhdx.support import np_from_txt
import pickle
import os
from numpy.lib.recfunctions import stack_arrays, append_fields

import os
import numpy as np
from pyhdx.support import np_from_txt, fmt_export
from pyhdx import PeptideMasterTable, KineticsFitting
import pickle

from bokeh.models import ColumnDataSource

directory = os.path.dirname(__file__)
np.random.seed(43)


fpath = os.path.join(directory, 'test_data', 'simulated_data.csv')
data = np_from_txt(fpath, delimiter=',')
data['end'] += 1
uptake = np.zeros_like(data, dtype=float)
data = append_fields(data, 'uptake', data=uptake, usemask=False)
sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
start, end = 5, 45  # total span of protein (inc, inc)
nc_start, nc_end = 31, 34  # span of no coverage area (inc, inc)

pmt = PeptideMasterTable(data, drop_first=1, ignore_prolines=True, remove_nan=False)

uptake_corrected = pmt.data['scores'] * pmt.data['ex_residues']
pmt.data = append_fields(pmt.data, ['uptake_corrected'], [uptake_corrected])


states = pmt.groupby_state()
series = states['state1']
series.make_uniform()

kf = KineticsFitting(series, bounds=(1e-2, 800))

ctrl = Controller('template', ['asdf'])
ctrl.series = series

print(ctrl.series.tf_cov)

pickle_names = ['fit_simulated_wt_avg', 'fit_simulated_blocks', 'fit_simulated_rates', 'fit_simulated_pfact']
fit_results = []

for name in pickle_names:
    with open(os.path.join(directory, 'test_data', name + '.pick'), 'rb') as f:
        fr = pickle.load(f)
        fit_results.append(fr)

#
# fr1 = {'rates': fit_results[0].output, 'fitresult': fit_results[0]}
# ctrl.fit_results['fit1'] = fr1
# fr2 = {'rates': fit_results[1].output, 'fitresult': fit_results[1]}

i = 0
output = fit_results[i].output
dic = {name: output[name] for name in output.dtype.names}
dic['y'] = output['rate']
dic['color'] = np.full_like(output, fill_value='red', dtype='<U7')
cds = ColumnDataSource(dic)
ctrl.sources['fit1'] = cds
ctrl.fit_results['fit1'] = fit_results[i]
ctrl.param.trigger('sources')
ctrl.param.trigger('fit_results')


i = 1
output = fit_results[i].output
dic = {name: output[name] for name in output.dtype.names}
dic['y'] = output['rate']
dic['color'] = np.full_like(output, fill_value='blue', dtype='<U7')
cds = ColumnDataSource(dic)
ctrl.sources['fit2'] = cds
ctrl.fit_results['fit2'] = fit_results[i]
ctrl.param.trigger('sources')

ctrl.fit_control.param['do_fit2'].constant = False
ctrl.param.trigger('fit_results')

ctrl.tf_fit_control.epochs=10

ctrl.tf_fit_control._do_fitting()
ctrl.tf_fit_control.param['do_fit'].constant = False

ctrl.classification_panel.target = 'pfact'
ctrl.classification_panel._action_threshold()


# if __name__ == '__main__':
#     ctrl.serve()
