import datetime
import logging
from pathlib import Path

import numpy as np
import panel as pn
import torch

from pyhdx.config import cfg
from pyhdx.local_cluster import verify_cluster
from pyhdx.web.apps import main_app, rfu_app
from pyhdx.web.base import STATIC_DIR

APP_DICT = {"main": lambda: main_app()[1], "rfu": lambda: rfu_app()[1]}


def run_apps():
    np.random.seed(43)
    torch.manual_seed(43)

    scheduler_address = cfg.get("cluster", "scheduler_address")
    if not verify_cluster(scheduler_address):
        print(
            f"No valid Dask scheduler found at specified address: '{scheduler_address}'"
        )

    log_root_dir = Path.home() / ".pyhdx" / "logs"
    log_dir = log_root_dir / datetime.datetime.now().strftime("%Y%m%d")
    log_dir.mkdir(
        parents=True, exist_ok=True
    )  # catch error when log dir does not exist
    root_log = logging.getLogger("pyhdx")
    root_log.setLevel(logging.DEBUG)

    fh = logging.FileHandler(log_dir / "pyhdx_logs.txt")
    formatter = logging.Formatter(
        "%(asctime)s %(name)s [%(levelname)s]: %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    root_log.addHandler(fh)
    root_log.info("Starting PyHDX server")

    tornado_logger = logging.getLogger("tornado.application")
    fh = logging.FileHandler(log_dir / "tornado_logs.txt")
    formatter = logging.Formatter(
        "%(asctime)s %(name)s [%(levelname)s]: %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    fh.setFormatter(formatter)
    fh.setLevel(10)
    tornado_logger.addHandler(fh)

    # TODO Clean assets dir from pdb files
    Path(cfg.assets_dir).mkdir(exist_ok=True, parents=True)

    print("Welcome to the PyHDX server!")
    pn.serve(
        APP_DICT,
        static_dirs={"pyhdx": STATIC_DIR, "assets": str(cfg.assets_dir)},
        index=str(STATIC_DIR / "index.html"),
    )


if __name__ == "__main__":
    run_apps()
