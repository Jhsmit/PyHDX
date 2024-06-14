from pyhdx.__version__ import __version__

VERSION_STRING = f"PyHDX {__version__}"

try:
    from pyhdx.output import FitReport
except ModuleNotFoundError:
    pass
