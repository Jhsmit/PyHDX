import pathlib

from bokeh.core.properties import Int, String
from bokeh.layouts import column
from bokeh.models import HTMLBox

CUSTOM_TS = pathlib.Path(__file__).parent / "html_button_model.ts"
CUSTOM_TS_STR = str(CUSTOM_TS.resolve())

DEFAULT_OBJECT = "<button style='width:100%'>Click Me</button>"


class HTMLButton(HTMLBox):
    """Example implementation of a Custom Bokeh Model"""

    __implementation__ = CUSTOM_TS_STR

    object = String(default=DEFAULT_OBJECT)
    clicks = Int(default=0)
