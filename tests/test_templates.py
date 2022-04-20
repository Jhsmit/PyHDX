import numpy as np
import matplotlib as mpl
from pyhdx.support import rgb_to_hex
import importlib
from pathlib import Path


class TestTemplates(object):
    def test_templates_converions(self):
        # import templates and run
        template_dir = Path(__file__).parent  # / something

        # todo fix
        scripts = [
            f.stem
            for f in template_dir.iterdir()
            if f.stem.startswith("fig") and f.suffix == ".py"
        ]
        for s in scripts:
            print(s)
            # importlib.import_module(s)
