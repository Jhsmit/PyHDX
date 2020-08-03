from dask.distributed import LocalCluster
import time
import os
from pathlib import Path

os.chdir(os.path.expanduser("~"))
print(os.getcwd())

print()
os.chdir(Path(__file__).parent.parent)

from pyhdx.math import solve_nnls


if __name__ == '__main__':


    cluster = LocalCluster(scheduler_port=52123, n_workers=10)
    print(cluster)

    try:
        loop = True
        while loop:
            try:
                time.sleep(2)
                print(cluster)
            except KeyboardInterrupt:
                print('interrupted')
                loop = False
    finally:
        print('closing')
        cluster.close()
