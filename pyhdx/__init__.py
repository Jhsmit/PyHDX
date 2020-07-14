# -*- coding: utf-8 -*-

"""Top-level package for PyHDX."""

from .pyhdx import PeptideMasterTable, PeptideMeasurements, KineticsSeries, Coverage
from .fitting import KineticsFitting
from .fileIO import read_dynamx
from .output import Output

from pbr.version import VersionInfo
from pbr import git
from pkg_resources import DistributionNotFound
import os
from pathlib import Path

package_name = 'pyhdx'
try:
    info = VersionInfo(package_name)
    __version__ = info.version_string()
    __dev_version__ = info.version_string_with_vcs()

except:
    __version__ = '0.0.0'
    __dev_version__ = '0.0.0'
    print('Warning: Version number not found.')

__git_sha__ = git.get_git_short_sha()
git_str = f' ({__git_sha__})' if __git_sha__ else ''
VERSION_STRING = f'pyHDX version {__version__}, development version {__dev_version__}' + git_str
