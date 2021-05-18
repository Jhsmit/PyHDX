from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, KineticsSeries
from pyhdx.fitting import BatchFitting, KineticsFitting, fit_rates_weighted_average, fit_rates
from pyhdx.fileIO import csv_to_protein
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
    sequences = np.random.choice(np.unique(data['sequence']), size=10, replace=False)
    bools = np.isin(data['sequence'], sequences)
    selected_data = data[bools]
    st = KineticsSeries(selected_data)
    data_objs.append(st)

if __name__ == '__main__':
    client = Client('tcp://127.0.0.1:52123')

    # results = []
    # for d in data_objs:
    #     result = fit_rates_weighted_average(d)
    #     results.append(result)
    # results = compute(*results)
    # print(results)

    async def do_fitting(data_objs):
        # client = await Client('tcp://127.0.0.1:52123')
        # futures = client.map(fit_rates_weighted_average, data_objs)
        # results = client.gather(futures)

        results = []
        for d in data_objs:
            result = fit_rates_weighted_average(d)
            results.append(result)

        compute(results)

        return results

#    print(futures)
    result = asyncio.run(do_fitting(data_objs))
    print(result)

    #results = client.gather(futures)

    #print(results)