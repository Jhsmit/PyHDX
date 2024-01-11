"""Load a fit result from a directory directly ino the web interface"""
from pyhdx.fileIO import load_fitresult
from pathlib import Path
from pyhdx.web.apps import main_app
from pyhdx.web.base import STATIC_DIR
import panel as pn

current_dir = Path(__file__).parent
fit_result = load_fitresult(
    current_dir.parent / "tests" / "test_data" / "output" / "ecsecb_tetramer_dimer"
)

ctrl, tmpl = main_app()

src = ctrl.sources["main"]
for hdxm in fit_result.hdxm_set:
    src.add(hdxm, hdxm.name)
src.add(fit_result, "fitresult_1")

pn.serve(tmpl, show=True, static_dirs={"pyhdx": STATIC_DIR})
