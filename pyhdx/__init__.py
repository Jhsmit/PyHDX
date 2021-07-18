from .models import PeptideMasterTable, PeptideMeasurements, HDXMeasurement, Coverage
from .fileIO import read_dynamx
from pathlib import Path

from ._version import get_versions


__version__ = get_versions()['version']

VERSION_STRING = f'PyHDX {__version__}'


del get_versions
