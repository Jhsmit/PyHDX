from __future__ import annotations

import argparse
import asyncio
import time
from asyncio import Future
from typing import Callable, Iterable, Any

from dask.distributed import LocalCluster, Client
from distributed import connect

from pyhdx.config import cfg
from pyhdx.support import select_config


class DummyClient(object):
    """Object to use as dask Client-like object for doing local operations with
    the dask Client API.
    """

    @staticmethod
    def submit(func: Callable, *args: Any, **kwargs) -> Future:
        future = Future()
        future.set_result(func(*args))
        return future

    @staticmethod
    def map(func: Callable, *iterables: Iterable, **kwargs) -> list[Future]:
        futures = []
        for items in zip(*iterables):
            result = func(*items)
            future = Future()
            future.set_result(result)
            futures.append(future)

        return futures

    @staticmethod
    def gather(futures) -> list[Any]:
        return [future.result() for future in futures]


def default_client(timeout="2s", **kwargs):
    """Return Dask client at scheduler adress as defined by the global config"""
    scheduler_address = cfg.cluster.scheduler_address
    try:
        client = Client(scheduler_address, timeout=timeout, **kwargs)
        return client
    except (TimeoutError, IOError):
        print(f"No valid Dask scheduler found at specified address: '{scheduler_address}'")
        return False


def default_cluster(**kwargs):
    """Start a dask LocalCluster at the scheduler port given by the config

    kwargs: override defaults

    """

    scheduler_address = cfg.cluster.scheduler_address
    port = int(scheduler_address.split(":")[-1])

    settings = {
        "scheduler_port": port,
        "n_workers": cfg.cluster.n_workers,
    }
    settings.update(kwargs)
    cluster = LocalCluster(**settings)

    return cluster


def verify_cluster(scheduler_address, timeout="2s"):
    """Check if a valid dask scheduler is running at the provided scheduler_address"""
    try:
        asyncio.run(connect(scheduler_address, timeout=timeout))
        return True
    except (TimeoutError, OSError):
        return False


def verify_cluster_async(scheduler_address, timeout="2s"):
    """Check if a valid dask scheduler is running at the provided scheduler_address"""
    try:
        asyncio.run(connect(scheduler_address, timeout=timeout))
        return True
    except (TimeoutError, OSError):
        return False


def blocking_cluster():
    """Start a dask LocalCluster and block until iterrupted"""
    parser = argparse.ArgumentParser(description="Start a new Dask local cluster")
    parser.add_argument("-p", "--port", help="Port to use for the Dask local cluster", dest="port")

    args = parser.parse_args()

    if args.port:
        port = int(args.port)
    else:
        scheduler_address = cfg.cluster.scheduler_address
        port = int(scheduler_address.split(":")[-1])
    try:
        n_workers = cfg.cluster.n_workers
        local_cluster = LocalCluster(scheduler_port=port, n_workers=n_workers)
        print(f"Started local cluster at {local_cluster.scheduler_address}")
    except OSError:
        print(f"Could not start local cluster with at port: {port}")
        raise
    try:
        loop = True
        while loop:
            try:
                time.sleep(2)
            except KeyboardInterrupt:
                print("Interrupted")
                loop = False
    finally:
        local_cluster.close()


if __name__ == "__main__":
    select_config()
    blocking_cluster()
