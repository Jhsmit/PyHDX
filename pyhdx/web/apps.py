
import yaml
from pathlib import Path
import collections

from pyhdx import VERSION_STRING
from pyhdx.batch_processing import yaml_to_hdxm
from pyhdx.web.base import STATIC_DIR
from pyhdx.web.constructor import AppConstructor, datetime
from pyhdx.fitting import RatesFitResult, GenericFitResult
from pyhdx.web.log import get_default_handler, StreamToLogger
from pyhdx.fileIO import csv_to_dataframe

import panel as pn
import sys

from pyhdx.web.log import logger

@logger('pyhdx')
def main_app():

    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load((cwd / 'pyhdx_app.yaml').read_text(encoding='utf-8'))

    ctr = AppConstructor(logger=main_app.logger)
    self = ctr

    ctr.parse(yaml_dict)
    ctrl = ctr.make_ctrl()

    ctrl.start()

    fmt = {'accent_base_color': '#1d417a'}

    tmpl = pn.template.FastGridTemplate(title=f'{VERSION_STRING}', **fmt)
    controllers = ctrl.control_panels.values()
    controls = pn.Accordion(*[controller.panel for controller in controllers], toggle=True)
    tmpl.sidebar.append(controls)

    views_names = [
        'rfu_scatter',
        'coverage',
        'logging_info',
        'logging_debug',
        'protein',
        'ddG',
        'rates',
        'gibbs'
        ]

    views = {v: ctrl.views[v] for v in views_names}
    [v.update() for v in views.values()]

    cov_tab = pn.Tabs(
        ('Coverage', views['coverage'].panel),
        ('Protein', views['protein'].panel)
    )

    scatter_tab = pn.Tabs(
        ('RFU', views['rfu_scatter'].panel),
        ('Rates', views['rates'].panel),
        ('ΔG', views['gibbs'].panel),
        ('ΔΔG', views['ddG'].panel),
    )

    log_tab = pn.Tabs(
        ('Info log', views['logging_info'].panel),
        ('Debug log', views['logging_debug'].panel)
    )

    tmpl.main[0:3, 0:6] = cov_tab
    tmpl.main[0:3, 6:12] = scatter_tab
    tmpl.main[3:5, 0:6] = log_tab

    return ctrl, tmpl