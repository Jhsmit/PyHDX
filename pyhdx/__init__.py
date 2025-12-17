from pyhdx.__version__ import __version__
from pyhdx.fitting_torch import TorchFitResult, TorchFitResultSet
from pyhdx.models import (
    Coverage,
    HDXMeasurement,
    HDXMeasurementSet,
    HDXTimepoint,
)

VERSION_STRING = f"PyHDX {__version__}"

try:
    from pyhdx.output import FitReport
except ModuleNotFoundError:
    pass


__all__ = [
    "HDXTimepoint",
    "HDXMeasurement",
    "Coverage",
    "HDXMeasurementSet",
    "TorchFitResult",
    "TorchFitResultSet",
    "FitReport",
]
