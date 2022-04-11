from pathlib import Path

import panel as pn
import yaml

from pyhdx import VERSION_STRING
from pyhdx.web.constructor import AppConstructor
from pyhdx.web.log import logger
from pyhdx.web.cache import MemoryCache, HybridHDFCache
from pyhdx.web.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.web.theme import ExtendedGoldenDefaultTheme, ExtendedGoldenDarkTheme

cache = MemoryCache(max_items=2000)


fmt = {"header_background": "#1d417a", "header_color": "#1d417a"}

# default kwargs for normal views (not the logs)
view_kwargs = {"scrollable": False}

fmt_kwargs = {**fmt}


@logger("pyhdx")
def main_app():
    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load(
        (cwd / "apps" / "pyhdx_app.yaml").read_text(encoding="utf-8")
    )

    ctr = AppConstructor(loggers={"pyhdx": main_app.logger}, cache=cache)

    ctrl = ctr.parse(yaml_dict)

    elvis = GoldenElvis(
        ctrl, ExtendedGoldenTemplate, ExtendedGoldenDefaultTheme, title=VERSION_STRING
    )
    tmpl = elvis.compose(
        elvis.column(
            elvis.row(  # top row
                elvis.stack(
                    elvis.view("coverage", **view_kwargs),
                    elvis.view("protein", **view_kwargs),
                    elvis.view("peptide_mse", title="Peptide MSE", **view_kwargs),
                    width=61.803,
                ),
                elvis.stack(
                    elvis.view("rfu_scatter", title="RFU", **view_kwargs),
                    elvis.view("drfu", title="ΔRFU", **view_kwargs),
                    elvis.view("rates", **view_kwargs),
                    elvis.view("gibbs_overlay", title="ΔG", **view_kwargs),
                    elvis.view("ddG_overlay", title="ΔΔG", **view_kwargs),
                ),
                height=61.803,
            ),
            elvis.row(  # second row
                elvis.stack(
                    elvis.view("logging_info", title="Info Log"),
                    elvis.view("logging_debug", title="Debug Log"),
                ),
                elvis.stack(
                    elvis.view("peptide_overlay", title="Peptide", **view_kwargs),
                    elvis.view("loss_lines", title="Losses", **view_kwargs),
                ),
            ),
        ),
        **fmt_kwargs
    )

    return ctrl, tmpl


@logger("pyhdx")
def rfu_app():
    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load(
        (cwd / "apps" / "rfu_app.yaml").read_text(encoding="utf-8")
    )

    ctr = AppConstructor(loggers={"pyhdx": rfu_app.logger}, cache=cache)
    ctrl = ctr.parse(yaml_dict)

    elvis = GoldenElvis(
        ctrl, ExtendedGoldenTemplate, ExtendedGoldenDefaultTheme, title=VERSION_STRING
    )

    tmpl = elvis.compose(
        elvis.column(
            elvis.row(  # top row
                elvis.stack(
                    elvis.view("coverage", title="Coverage", **view_kwargs),
                    elvis.view("protein", **view_kwargs),
                ),
                elvis.stack(
                    elvis.view("rfu_scatter", title="RFU", **view_kwargs),
                    elvis.view("drfu", title="ΔRFU", **view_kwargs),
                ),
            ),
            elvis.row(  # second row
                elvis.stack(
                    elvis.view("logging_info", title="Info Log"),
                    elvis.view("logging_debug", title="Debug Log"),
                ),
                elvis.stack(
                    elvis.view("peptide_scatter", title="Peptide", **view_kwargs),
                ),
            ),
        ),
        **fmt_kwargs
    )

    return ctrl, tmpl
