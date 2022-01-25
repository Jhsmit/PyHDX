import pathlib

import panel as pn
import param

from pyhdx.web.template import ExtendedGoldenTemplate

from panel.template.fast.theme import DEFAULT_STYLE, DARK_STYLE
from bokeh.themes import Theme as BkTheme

custom_default_json_theme = DEFAULT_STYLE.create_bokeh_theme()
custom_default_json_theme["attrs"]["ColorBar"]["background_fill_alpha"] = 0
custom_dark_json_theme = DARK_STYLE.create_bokeh_theme()
custom_dark_json_theme["attrs"]["ColorBar"]["background_fill_alpha"] = 0


class ExtendedGoldenDefaultTheme(pn.template.golden.GoldenDefaultTheme):

    css = param.Filename(
        default=pathlib.Path(__file__).parent
        / "static"
        / "extendedgoldentemplate"
        / "default.css"
    )

    _template = ExtendedGoldenTemplate

    bokeh_theme = param.ClassSelector(
        class_=(BkTheme, str), default=BkTheme(json=custom_default_json_theme)
    )


class ExtendedGoldenDarkTheme(pn.template.golden.GoldenDarkTheme):

    css = param.Filename(
        default=pathlib.Path(__file__).parent
        / "static"
        / "extendedgoldentemplate"
        / "dark.css"
    )

    _template = ExtendedGoldenTemplate

    bokeh_theme = param.ClassSelector(
        class_=(BkTheme, str), default=BkTheme(json=custom_dark_json_theme)
    )
