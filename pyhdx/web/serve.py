import panel as pn
from pyhdx.web.apps import main_app#, diff_app, single_app, folding_app, full_deuteration_app
from pyhdx.web.base import STATIC_DIR
import numpy as np
import torch

from pyhdx.config import ConfigurationSettings
from pyhdx.local_cluster import verify_cluster

import logging
from pathlib import Path
import datetime

APP_DICT = {
    'main': lambda: main_app().template,
    # 'single': lambda: single_app().template,
    # 'diff': lambda: diff_app().template,
    # 'folding': lambda: folding_app().template,
    # 'full_deuteration': lambda: full_deuteration_app().template
}


def run_main():
    np.random.seed(43)
    torch.manual_seed(43)

    scheduler_address = ConfigurationSettings().get('cluster', 'scheduler_address')
    if not verify_cluster(scheduler_address):
        print(f"No valid Dask scheduler found at specified address: '{scheduler_address}'")
        return

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
    pn.serve(APP_DICT, static_dirs={'pyhdx': STATIC_DIR})


if __name__ == '__main__':
    run_main()


