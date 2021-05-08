from dask.distributed import LocalCluster
import time
import os
from pathlib import Path
from pyhdx.panel.config import ConfigurationSettings
import argparse

cfg = ConfigurationSettings()


def default_cluster():
    port = int(cfg.get('cluster', 'port'))
    cluster = LocalCluster(scheduler_port=port, n_workers=10)  # todo default settings local cluster from config

    return cluster


def blocking_cluster():
    parser = argparse.ArgumentParser(description='Start a new Dask local cluster')
    parser.add_argument('-p', '--port', help="Port to use for the Dask local cluster", dest='port')

    args = parser.parse_args()

    if args.port:
        port = int(args.port)
    else:
        port = int(cfg.get('cluster', 'port'))

    local_cluster = LocalCluster(scheduler_port=port, n_workers=10)  # todo default settings local cluster from config
    print(f"Started local cluster at {local_cluster.scheduler_address}")
    try:
        loop = True
        while loop:
            try:
                time.sleep(2)
            except KeyboardInterrupt:
                print('Interrupted')
                loop = False
    finally:
        local_cluster.close()


if __name__ == '__main__':
    # import sys
    # sys.argv.append('-p')
    # sys.argv.append('52348')
    blocking_cluster()
