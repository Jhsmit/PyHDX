import panel as pn
from pyhdx.panel.apps import main_app, diff_app, single_app, folding_app, full_deuteration_app
from pyhdx.panel.base import STATIC_DIR
import numpy as np
import torch

APP_DICT = {
    'main': main_app,
    'single': single_app,
    'diff': diff_app,
    'folding': folding_app,
    'full_deuteration': full_deuteration_app
}

def run_main():
    np.random.seed(43)
    torch.manual_seed(43)
    pn.serve(APP_DICT, static_dirs={'pyhdx': STATIC_DIR})

if __name__ == '__main__':
    run_main()


