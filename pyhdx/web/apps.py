from pathlib import Path

import panel as pn
import param
import yaml
import re

from pyhdx import VERSION_STRING
from pyhdx.web.constructor import AppConstructor
from pyhdx.web.log import logger
from pyhdx.web.cache import MemoryCache
from pyhdx.web.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.web.theme import ExtendedGoldenDefaultTheme

cache = MemoryCache(max_items=2000)

# Check for new panel releases if this is still needed
pn.extension("mathjax")

fmt = {"header_background": "#1d417a", "header_color": "#1d417a"}

# default kwargs for normal views (not the logs)
view_kwargs = {"scrollable": False}

fmt_kwargs = {**fmt}

yaml.SafeLoader.add_constructor("!regexp", lambda l, n: re.compile(l.construct_scalar(n)))


@logger("pyhdx")
def main_app():
    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load((cwd / "apps" / "pyhdx_app.yaml").read_text(encoding="utf-8"))

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
                    width=50,
                ),
                elvis.stack(
                    elvis.view("rfu_overlay", title="RFU", **view_kwargs),
                    elvis.view("drfu_overlay", title="ΔRFU", **view_kwargs),
                    elvis.view("d_uptake_overlay", title="D-Uptake", **view_kwargs),
                    elvis.view("dd_uptake_scatter", title="ΔD-Uptake", **view_kwargs),
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
    yaml_dict = yaml.safe_load((cwd / "apps" / "rfu_app.yaml").read_text(encoding="utf-8"))

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
                    elvis.view("rfu_overlay", title="RFU", **view_kwargs),
                    elvis.view("drfu_overlay", title="ΔRFU", **view_kwargs),
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


# These classes should be moved elsewhere
class SliderLayout(param.Parameterized):
    def __init__(self, ctrl, **params):
        super().__init__(**params)
        self.ctrl = ctrl

        self.dG_text = pn.widgets.StaticText(value="<b>ΔG</b> (kJ/mol)")
        self.k_open_text = pn.widgets.StaticText(value="<b>k_open</b> (Log10, s⁻¹)")
        self.k_close_text = pn.widgets.StaticText(value="<b>k_close</b> (Log10, s⁻¹)")
        self.layout = pn.Column()
        self.update()

        self.ctrl.param.watch(self.update, ["updated"])

    def update(self, *events) -> None:
        self.layout[:] = [
            self.dG_text,
            self.ctrl.widgets["dG"],
            self.k_open_text,
            self.ctrl.widgets["k_open"],
            self.k_close_text,
            self.ctrl.widgets["k_close"],
        ]

    @property
    def panel(self) -> pn.Column:
        return self.layout


class MarkDownFile(param.Parameterized):
    def __init__(self, pth, **params):
        super().__init__(**params)
        self.object = Path(pth).read_text(encoding="utf-8")

    @property
    def panel(self):
        style = {"font-size": "16pt"}

        return pn.pane.Markdown(
            self.object,
            sizing_mode="stretch_both",
            disable_math=False,
            style=style,
        )


@logger("pyhdx")
def peptide_app():
    cwd = Path(__file__).parent.resolve()
    yaml_dict = yaml.safe_load((cwd / "apps" / "peptide_app.yaml").read_text(encoding="utf-8"))

    ctr = AppConstructor(cache=cache)

    ctrl = ctr.parse(yaml_dict)
    peptide_ctrl = ctrl.control_panels["PeptidePropertiesControl"]

    elvis = GoldenElvis(
        ctrl, ExtendedGoldenTemplate, ExtendedGoldenDefaultTheme, title=VERSION_STRING
    )

    slider_layout = SliderLayout(peptide_ctrl)

    faq_path = cwd / "apps" / "peptide_faq.md"
    md_pane = MarkDownFile(faq_path)

    tmpl = elvis.compose(
        elvis.row(
            elvis.column(  # left column
                elvis.view("peptide_uptake", title="Peptide D-uptake", **view_kwargs),
                elvis.stack(
                    elvis.view("aa_uptake", title="AA D-uptake", **view_kwargs),
                    elvis.view("k_int", title="Intrinsic rates", **view_kwargs),
                ),
                width=61.803,
            ),
            elvis.column(  # right column
                elvis.stack(
                    elvis.view(slider_layout, title="Rate sliders"),
                    elvis.view(md_pane, title="FAQ"),
                )
            ),
        ),
        **fmt_kwargs
    )

    return ctrl, tmpl
