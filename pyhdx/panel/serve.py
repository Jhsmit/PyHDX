import panel as pn
from pyhdx.panel.apps import main_app#, diff_app, single_app, folding_app, full_deuteration_app
from pyhdx.panel.base import STATIC_DIR
import numpy as np
import torch

from pyhdx.panel.config import ConfigurationSettings
from pyhdx.support import verify_cluster

import logging
from functools import partial
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

    cluster = ConfigurationSettings().cluster
    if verify_cluster(cluster):
        print("Welcome to the PyHDX server!")
        pn.serve(APP_DICT, static_dirs={'pyhdx': STATIC_DIR})
    #todo add some logic/kwargs for log keeping
    log_root_dir = Path('logs')
    log_dir = log_root_dir / datetime.datetime.now().strftime('%Y%m%d')
    log_dir.mkdir(parents=False, exist_ok=True) # catch error when log dir does not exist

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

    #sys.stderr = StreamToLogger(logger, logging.DEBUG)

    pn.serve(APP_DICT, static_dirs={'pyhdx': STATIC_DIR})


if __name__ == '__main__':
    root_log = logging.getLogger('pyhdx')
    run_main()


