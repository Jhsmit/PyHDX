from .models import (
    HDXTimepoint,
    HDXMeasurement,
    Coverage,
    HDXMeasurementSet,
)
from .fileIO import read_dynamx
from .fitting_torch import TorchFitResult, TorchFitResultSet
from .batch_processing import StateParser
from ._version import get_versions

try:
    from .output import FitReport
except ModuleNotFoundError:
    pass


__version__ = get_versions()["version"]

VERSION_STRING = f"PyHDX {__version__}"


del get_versions
