# -*- coding: utf-8 -*-

"""Top-level package for PyHDX."""

from .pyhdx import PeptideMasterTable, PeptideMeasurements, KineticsSeries, Coverage
from .fitting import KineticsFitting
from .fileIO import read_dynamx
from .output import Output


from pbr.version import VersionInfo
package_name = 'pyhdx'
info = VersionInfo(package_name)
__version__ = info.version_string()
