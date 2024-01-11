from pathlib import Path
import sys


class TestTemplates(object):
    def test_templates_converions(self):
        # import templates and run
        template_dir = Path(__file__).parent.parent / "templates"  # / something
        print(template_dir)
        sys.path.append(str(template_dir))

        scripts = [
            f.stem for f in template_dir.iterdir() if f.stem[:2].isdigit() and f.suffix == ".py"
        ]
        for s in scripts:
            pass
            # doesnt work with some scripts requiring a cluster
            # print(s)
            # importlib.import_module(s)
