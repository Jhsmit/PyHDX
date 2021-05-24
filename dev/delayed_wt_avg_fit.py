from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, KineticsSeries
from pyhdx.fitting import BatchFitting, KineticsFitting, fit_rates_weighted_average, fit_rates
from pyhdx.fileIO import csv_to_protein
from pyhdx.local_cluster import default_client
from dask.distributed import Client
from dask import delayed, compute
import numpy as np
import asyncio

current_dir = Path(__file__).parent
np.random.seed(43)

data_dir = current_dir.parent / 'tests' / 'test_data'
data = read_dynamx(data_dir / 'ecSecB_apo.csv', data_dir / 'ecSecB_dimer.csv')

pmt = PeptideMasterTable(data)
pmt.set_control(('Full deuteration control', 0.167))

names = ['SecB his dimer apo', 'SecB WT apo']
data_objs = []
for name in names:
    data = pmt.get_state(name)
    bools = data['end'] < 40

    selected_data = data[bools]
    st = KineticsSeries(selected_data)
    #st = KineticsSeries(data)

    data_objs.append(st)


if __name__ == '__main__':
    client = default_client()

    futures = client.map(fit_rates_weighted_average, data_objs, client='worker_client')
    print(futures)
    def future_done(future):
        print('joohoehoeho', future)
    futures[0].add_done_callback(future_done)

    result = client.gather(futures)
    print(result)

    print(result[0].output)