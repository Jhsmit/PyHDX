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

directory = os.path.dirname(__file__)
np.random.seed(43)


fpath = os.path.join(directory, 'test_data', 'simulated_data.csv')
data = np_from_txt(fpath, delimiter=',')
uptake = np.zeros_like(data, dtype=float)
data = append_fields(data, 'uptake', data=uptake, usemask=False)
sequence = 'XXXXTPPRILALSAPLTTMMFSASALAPKIXXXXLVIPWINGDKG'

timepoints = [0.167, 0.5, 1, 5, 10, 30, 100]
start, end = 5, 45  # total span of protein (inc, inc)
nc_start, nc_end = 31, 34  # span of no coverage area (inc, inc)

pmt = PeptideMasterTable(data, drop_first=0, ignore_prolines=False, remove_nan=False)
states = pmt.groupby_state()
series = states['state1']

kf = KineticsFitting(series, bounds=(1e-2, 800))

ctrl = Controller('template', ['asdf'])
ctrl.series = series

pickle_names = ['fit_simulated_wt_avg', 'fit_simulated_blocks', 'fit_simulated_rates', 'fit_simulated_pfact']
fit_results = []

for name in pickle_names:
    with open(os.path.join(directory, 'test_data', name + '.pick'), 'rb') as f:
        fr = pickle.load(f)
        fit_results.append(fr)

#
fr1 = {'rates': fit_results[0].output, 'fitresult': fit_results[0]}
ctrl.fit_results['fit1'] = fr1
fr2 = {'rates': fit_results[1].output, 'fitresult': fit_results[1]}

ctrl.fit_results['fit2'] = fr2
ctrl.fit_control.param['do_fit2'].constant = False
ctrl.param.trigger('fit_results')


#
# ctrl.tf_fit_control.param['do_fit'].constant = False

if __name__ == '__main__':
    ctrl.serve()
