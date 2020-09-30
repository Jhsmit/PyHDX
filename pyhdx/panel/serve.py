import panel as pn
from pyhdx.panel.apps import main_app, diff_app, single_app

APP_DICT = {
    'main': main_app,
    'single': diff_app,
    'diff': single_app
}


pn.serve(APP_DICT)


