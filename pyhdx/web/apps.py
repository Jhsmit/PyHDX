
from pathlib import Path

import panel as pn
import yaml

from pyhdx import VERSION_STRING
from pyhdx.web.constructor import AppConstructor
from pyhdx.web.log import logger


@logger('pyhdx')
def main_app():
    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load((cwd / 'pyhdx_app.yaml').read_text(encoding='utf-8'))

    ctr = AppConstructor(loggers={'pyhdx': main_app.logger})

    ctrl = ctr.parse(yaml_dict)

    ctrl.start()

    fmt = {
        'header_color': '#ffffff',  # this is the text
        'header_background': '#00407A',
        'accent_base_color': '#00407A',
        'theme_toggle': False
    }

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
        'ddG_overlay',
        'rates',
        'gibbs_overlay',
        'peptide_mse',
        #'peptide_scatter',
        'peptide_overlay',
        'loss_lines'
        ]

    views = {v: ctrl.views[v] for v in views_names}
    [v.update() for v in views.values()]

    # this should be on the view intances probably
    def get_view(name):
        return pn.Column(views[name].panel, sizing_mode='stretch_both')

    cov_tab = pn.Tabs(
        ('Coverage', get_view('coverage')),
        ('Protein', get_view('protein')),
        ('Peptide MSE', get_view('peptide_mse'))
    )

    scatter_tab = pn.Tabs(
        ('RFU', get_view('rfu_scatter')),
        ('Rates', get_view('rates')),
        ('ΔG', get_view('gibbs_overlay')),
        ('ΔΔG', get_view('ddG_overlay')),
    )

    log_tab = pn.Tabs(
        ('Info log', views['logging_info'].panel),
        ('Debug log', views['logging_debug'].panel)
    )

    peptide_tab = pn.Tabs(
        ('Peptide', get_view('peptide_overlay')),
        ('Losses', get_view('loss_lines'))
    )

    tmpl.main[0:3, 0:6] = cov_tab
    tmpl.main[0:3, 6:12] = scatter_tab
    tmpl.main[3:5, 0:6] = log_tab
    tmpl.main[3:5, 6:12] = peptide_tab

    return ctrl, tmpl