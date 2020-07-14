# -*- coding: utf-8 -*-

"""Top-level package for PyHDX."""

from .pyhdx import PeptideMasterTable, PeptideMeasurements, KineticsSeries, Coverage
from .fitting import KineticsFitting
from .fileIO import read_dynamx
from .output import Output

from pbr.version import VersionInfo
from pbr import git

package_name = 'pyhdx'
info = VersionInfo(package_name)
__version__ = info.version_string()
__dev_version__ = info.version_string_with_vcs()
__git_sha__ = git.get_git_short_sha()

