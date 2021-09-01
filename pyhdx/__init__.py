from .models import PeptideMasterTable, PeptideMeasurements, HDXMeasurement, Coverage, HDXMeasurementSet
from .fileIO import read_dynamx
from .fitting_torch import TorchSingleFitResult, TorchBatchFitResult
from ._version import get_versions


__version__ = get_versions()['version']

VERSION_STRING = f'PyHDX {__version__}'


del get_versions
