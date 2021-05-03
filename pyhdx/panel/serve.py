import panel as pn
from pyhdx.panel.apps import main_app, diff_app, single_app, folding_app, full_deuteration_app
from pyhdx.panel.base import STATIC_DIR
import numpy as np
import torch
from pyhdx.panel.config import ConfigurationSettings
from pyhdx.support import verify_cluster

APP_DICT = {
    'main': lambda: main_app().template,
    'single': lambda: single_app().template,
    'diff': lambda: diff_app().template,
    'folding': lambda: folding_app().template,
    'full_deuteration': lambda: full_deuteration_app().template
}


def run_main():
    np.random.seed(43)
    torch.manual_seed(43)

    cluster = ConfigurationSettings().cluster
    if verify_cluster(cluster):
        print("Welcome to the PyHDX server!")
        pn.serve(APP_DICT, static_dirs={'pyhdx': STATIC_DIR})


if __name__ == '__main__':
    run_main()


