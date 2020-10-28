import param
import panel as pn
import pathlib

from pyhdx.panel.template import ExtendedGoldenTemplate


class ExtendedGoldenDefaultTheme(pn.template.golden.GoldenDefaultTheme):

    css = param.Filename(default=pathlib.Path(__file__).parent / 'static' / 'extendedgoldentemplate' / 'default.css')

    _template = ExtendedGoldenTemplate


class ExtendedGoldenDarkTheme(pn.template.golden.GoldenDarkTheme):

    css = param.Filename(default=pathlib.Path(__file__).parent / 'static' / 'extendedgoldentemplate' / 'dark.css')

    _template = ExtendedGoldenTemplate
