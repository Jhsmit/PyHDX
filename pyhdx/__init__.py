"""Top-level package for PyHDX."""

from .models import PeptideMasterTable, PeptideMeasurements, KineticsSeries, Coverage
from .fitting import KineticsFitting
from .fileIO import read_dynamx
from .output import Output

from pbr import version
from pbr import git
from pkg_resources import DistributionNotFound
import os
from pathlib import Path


package_name = 'pyhdx'
try:
    info = version.VersionInfo(package_name)
    __version__ = info.version_string()
    __dev_version__ = info.version_string_with_vcs()
    __git_sha__ = git.get_git_short_sha()


except Exception:  # Pbr throws very broad Exception, for some reason DistributionNotFound does not want to be caught
    git_dir = Path(__file__).parent.parent / '.git'
    try:
        tagged = git._run_git_command(
            ['describe', '--tags'], git_dir,
            throw_on_error=True).replace('-', '.')
        semantic_version = version.SemanticVersion.from_pip_string(tagged)
        __version__ = semantic_version.brief_string()
        __dev_version__ = semantic_version._long_version(None)
        __git_sha__ = git.get_git_short_sha(git_dir)
    except FileNotFoundError:
        print('git not installed')
        __version__ = '0.0.0'
        __dev_version__ = '0.0.0'
        __git_sha__ = ''

git_str = f' ({__git_sha__})' if __git_sha__ else ''
VERSION_STRING = f'pyHDX version {__version__}, development version {__dev_version__}' + git_str
VERSION_STRING_SHORT = f'pyHDX v{__version__} ({__dev_version__}, {git_str})'
