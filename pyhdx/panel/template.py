import pathlib

from panel.template import GoldenTemplate

print(pathlib.Path(__file__).parent)

class ExtendedGoldenTemplate(GoldenTemplate):

    _template = pathlib.Path(__file__).parent / 'golden.html'

#    _css = pathlib.Path(__file__).parent / 'golden.css'
