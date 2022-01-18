from pathlib import Path

import panel as pn
import yaml

from pyhdx import VERSION_STRING
from pyhdx.web.constructor import AppConstructor
from pyhdx.web.log import logger
from pyhdx.web.cache import MemoryCache, HybridHDFCache
from pyhdx.web.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.web.theme import ExtendedGoldenDefaultTheme

cache = MemoryCache(max_items=2000)

#cache = HybridHDFCache(file_path ='test123.h5')


fmt = {
    'header_color': '#ffffff',  # this is the text
    'header_background': '#00407A',
    'accent_base_color': '#00407A',
    'theme_toggle': False
}

@logger('pyhdx')
def main_app():
    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load((cwd / 'apps' / 'pyhdx_app.yaml').read_text(encoding='utf-8'))

    ctr = AppConstructor(loggers={'pyhdx': main_app.logger}, cache=cache)

    ctrl = ctr.parse(yaml_dict)

    ctrl.start()

    tmpl = pn.template.FastGridTemplate(title=f'{VERSION_STRING}', **fmt)
    controllers = ctrl.control_panels.values()
    controls = pn.Accordion(*[controller.panel for controller in controllers], toggle=True)
    tmpl.sidebar.append(controls)

    elvis = GoldenElvis(ctrl, ExtendedGoldenTemplate, ExtendedGoldenDefaultTheme,
                        title=VERSION_STRING)

    tmpl = elvis.compose(
        elvis.column(
            elvis.row(  # top row
                elvis.stack(
                    elvis.view('coverage'),
                    elvis.view('protein'),
                    elvis.view('peptide_mse', title='Peptide MSE')
                ),
                elvis.stack(
                    elvis.view('rfu_scatter', title='RFU'),
                    elvis.view('drfu', title='ΔRFU'),
                    elvis.view('rates'),
                    elvis.view('gibbs_overlay', title='ΔG'),
                    elvis.view('ddG_overlay', title='ΔΔG')
                )
            ),
            elvis.row(  # second row
                elvis.stack(
                    elvis.view('logging_info', title='Info Log'),
                    elvis.view('logging_debug', title='Debug Log')
                ),
                elvis.stack(
                    elvis.view('peptide_overlay', title='Peptide'),
                    elvis.view('loss_lines', title='Losses')
                )
            )
        )
    )

    return ctrl, tmpl


@logger('pyhdx')
def rfu_app():
    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load((cwd / 'apps' / 'rfu_app.yaml').read_text(encoding='utf-8'))

    ctr = AppConstructor(loggers={'pyhdx': rfu_app.logger}, cache=cache)
    ctrl = ctr.parse(yaml_dict)

    ctrl.start()

    elvis = GoldenElvis(ctrl, ExtendedGoldenTemplate, ExtendedGoldenDefaultTheme,
                        title=VERSION_STRING)

    tmpl = elvis.compose(
        elvis.column(
            elvis.row(  # top row
                elvis.stack(
                    elvis.view('coverage', title='Coverage'),
                    elvis.view('protein'),
                ),
                elvis.stack(
                    elvis.view('rfu_scatter', title='RFU'),
                    elvis.view('drfu', title='ΔRFU'),
                )
            ),
            elvis.row(  # second row
                elvis.stack(
                    elvis.view('logging_info', title='Info Log'),
                    elvis.view('logging_debug', title='Debug Log')
                ),
                elvis.stack(
                    elvis.view('peptide_scatter', title='Peptide'),
                )
            )
        )
    )

    return ctrl, tmpl