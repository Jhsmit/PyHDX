import datetime
import logging
from pathlib import Path

import numpy as np
import panel as pn
import torch

from pyhdx.config import cfg
from pyhdx.local_cluster import verify_cluster, default_client
from pyhdx.web.apps import main_app, rfu_app
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.cache import MemoryCache, HybridHDFCache




def run_apps(executor, cache):
    np.random.seed(43)
    torch.manual_seed(43)

    APP_DICT = {  #TODO perhaps they should take a constructor instance?
        'main': lambda: main_app(executor=executor, cache=cache)[1],
        'rfu': lambda: rfu_app(executor=executor, cache=cache)[1]
    }


    log_root_dir = Path.home() / '.pyhdx' / 'logs'
    log_dir = log_root_dir / datetime.datetime.now().strftime('%Y%m%d')
    log_dir.mkdir(parents=True, exist_ok=True) # catch error when log dir does not exist
    root_log = logging.getLogger('pyhdx')
    root_log.setLevel(logging.DEBUG)

    fh = logging.FileHandler(log_dir / 'pyhdx_logs.txt')
    formatter = logging.Formatter('%(asctime)s %(name)s [%(levelname)s]: %(message)s', "%Y-%m-%d %H:%M:%S")
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    root_log.addHandler(fh)
    root_log.info('Starting PyHDX server')

    tornado_logger = logging.getLogger('tornado.application')
    fh = logging.FileHandler(log_dir / 'tornado_logs.txt')
    formatter = logging.Formatter('%(asctime)s %(name)s [%(levelname)s]: %(message)s', "%Y-%m-%d %H:%M:%S")
    fh.setFormatter(formatter)
    fh.setLevel(10)
    tornado_logger.addHandler(fh)

    print("Welcome to the PyHDX server!")
    pn.serve(APP_DICT, static_dirs={'pyhdx': STATIC_DIR}, index=str(STATIC_DIR / 'index.html'))


if __name__ == '__main__':

    scheduler_address = cfg.get('cluster', 'scheduler_address')
    if not verify_cluster(scheduler_address):
        print(f"No valid Dask scheduler found at specified address: '{scheduler_address}'")
    else:
        executor = default_client(asynchronous=True)
        cache = MemoryCache(max_items=2000)

        run_apps(executor, cache)


