import pathlib

from panel.template import GoldenTemplate


class ExtendedGoldenTemplate(GoldenTemplate):

    _template = pathlib.Path(__file__).parent / 'golden.html'

#    _css = pathlib.Path(__file__).parent / 'golden.css'
