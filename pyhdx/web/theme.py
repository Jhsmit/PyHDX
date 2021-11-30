import pathlib

import panel as pn
import param

from pyhdx.web.template import ExtendedGoldenTemplate


class ExtendedGoldenDefaultTheme(pn.template.golden.GoldenDefaultTheme):

    css = param.Filename(default=pathlib.Path(__file__).parent / 'static' / 'extendedgoldentemplate' / 'default.css')

    _template = ExtendedGoldenTemplate


class ExtendedGoldenDarkTheme(pn.template.golden.GoldenDarkTheme):

    css = param.Filename(default=pathlib.Path(__file__).parent / 'static' / 'extendedgoldentemplate' / 'dark.css')

    _template = ExtendedGoldenTemplate
