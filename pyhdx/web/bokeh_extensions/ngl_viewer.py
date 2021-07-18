import pathlib
from bokeh.core.properties import List, String, Bool
from bokeh.models import LayoutDOM


CUSTOM_TS = pathlib.Path(__file__).parent / "ngl_viewer.ts"
CUSTOM_TS_STR = str(CUSTOM_TS.resolve())


class ngl(LayoutDOM):
    __implementation__ = CUSTOM_TS_STR

    spin = Bool
    representation = String
    rcsb_id = String
    no_coverage = String
    color_list = List(List(String))
    pdb_string = String
