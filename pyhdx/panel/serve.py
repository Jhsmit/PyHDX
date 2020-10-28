import panel as pn
from pyhdx.panel.apps import main_app, diff_app, single_app, folding_app
from pyhdx.panel.base import STATIC_DIR

APP_DICT = {
    'main': main_app,
    'single': single_app,
    'diff': diff_app,
    'folding': folding_app
}


pn.serve(APP_DICT, static_dirs={'pyhdx': STATIC_DIR})


