from pyhdx import PeptideMeasurements, PeptideCSVFile, Coverage, KineticsSeries
from pyhdx.plot import plot_residue_map, make_coverage_figure
from pyhdx.fitting import KineticsFitting, KineticsFitResult, TwoComponentAssociationModel
import numpy as np
import matplotlib.pyplot as plt
#
#
filename = r"C:\Users\jhs\Programming\pyhdx\tests\test_data\ds1.csv"
control_100 = ('PpiA-FD', 0.167)
state = 'PpiANative'

pcf = PeptideCSVFile(filename)
# b = pcf.data['start'] < 50
# pcf.data = pcf.data[b]
print(len(pcf.data))
states_c = pcf.groupby_state_control(control_100)
series = states_c[state]

print('scf, scfA')
print(series[-1].scores)
print(series[-1].scores_average)


kf = KineticsFitting(series)

scores = kf.scores_norm
print(scores.shape)
#
# d = np.diff(scores, axis=1)
# bs = np.all(d == 0, axis=0)
# print(d)
# print(bs)
print(series[0].block_length)


m = BiexpIncreaseModel()

print(m)


m.initial_guess(series.times, scores[:, 0])

if __name__ == '__main__':
    fr = kf.do_fitting()

    res = fr.results[0]
    model = fr.models[0]
    print(res.model)
    print(fr.models[0])
    print(res.params)
    print('hoidoei')


    plt.figure()
    plt.plot(fr.rate)
    plt.yscale('log')
    plt.show()
