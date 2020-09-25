import panel as pn

from pyhdx.panel.main_app import main_app as main_ctrl
from pyhdx.panel.compare_app import compare_app as diff_ctrl
from pyhdx.panel.single_app import single_app as single_ctrl

APP_DICT = {
    'main': main_ctrl,
    'single': single_ctrl,
    'diff': diff_ctrl
}


port = 51337
pn.serve(APP_DICT, websocket_origin=[f'localhost:{port}', 'pyhdx.jhsmit.org'], port=port, threaded=True)


