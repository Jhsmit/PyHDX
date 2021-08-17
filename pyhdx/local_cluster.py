from dask.distributed import LocalCluster, Client
import time
from pyhdx.config import ConfigurationSettings
import argparse

cfg = ConfigurationSettings()

#todo refactor cluster
def default_client(timeout='2s'):
    scheduler_address = cfg.get('cluster', 'scheduler_address')
    try:
        client = Client(scheduler_address, timeout=timeout)
        return client
    except (TimeoutError, IOError):
        print(f"No valid Dask scheduler found at specified address: '{scheduler_address}'")
        return False


def default_cluster(**kwargs):
    """Start a dask LocalCluster at the scheduler port given by the config

    kwargs: override defaults

    """

    scheduler_address = cfg.get('cluster', 'scheduler_address')
    port = int(scheduler_address.split(':')[-1])

    settings = {
        'scheduler_port': port,
        'n_workers': int(cfg.get('cluster', 'n_workers'))}
    settings.update(kwargs)
    cluster = LocalCluster(**settings)

    return cluster


def verify_cluster(scheduler_address, timeout='2s'):
    """Check if a valid dask scheduler is running at the provided scheduler_address"""
    try:
        client = Client(scheduler_address, timeout=timeout)
        return True
    except (TimeoutError, IOError):
        return False


def blocking_cluster():
    """Start a dask LocalCluster and block until iterrupted"""
    parser = argparse.ArgumentParser(description='Start a new Dask local cluster')
    parser.add_argument('-p', '--port', help="Port to use for the Dask local cluster", dest='port')

    args = parser.parse_args()

    if args.port:
        port = int(args.port)
    else:
        scheduler_address = cfg.get('cluster', 'scheduler_address')
        port = int(scheduler_address.split(':')[-1])
    try:
        n_workers = int(cfg.get('cluster', 'n_workers'))
        local_cluster = LocalCluster(scheduler_port=port, n_workers=n_workers)
        print(f"Started local cluster at {local_cluster.scheduler_address}")
    except OSError as e:
        print(f"Could not start local cluster with at port: {port}")
        raise
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
