import panel as pn

from pyhdx.panel.main import tmpl as main_ctrl
from pyhdx.panel.apps import tmpl as diff_ctrl
from pyhdx.panel.single_app import tmpl as single_ctrl

APP_DICT = {
    'main': main_ctrl,
    'single': single_ctrl,
    'diff': diff_ctrl
}


port = 51337
pn.serve(APP_DICT, websocket_origin=[f'localhost:{port}', 'pyhdx.jhsmit.org'], port=port)


