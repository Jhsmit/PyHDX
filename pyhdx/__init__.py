"""Top-level package for PyHDX."""
import setuptools   # Import prevents warning
from .models import PeptideMasterTable, PeptideMeasurements, KineticsSeries, Coverage
from .fitting import KineticsFitting
from .fileIO import read_dynamx
from .output import Output
from pathlib import Path
import importlib.metadata


package_name = 'pyhdx'

__version__ = importlib.metadata.version(package_name)

try:
    from pbr import version
    from pbr import git
    has_pbr = True
except ModuleNotFoundError:
    has_pbr = False

info = version.VersionInfo(package_name)
git_dir = Path(__file__).parent.parent / '.git'

if has_pbr:
    try:
        __git_sha__ = git.get_git_short_sha(git_dir)

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
            # Git not installed
            __git_sha__ = None
else:
    __git_sha__ = None

VERSION_STRING = f'pyHDX version {__version__}'
VERSION_STRING_SHORT = f'pyHDX v{__version__}'
if __git_sha__ is not None:
    VERSION_STRING += f', development version {__git_sha__}'
    VERSION_STRING_SHORT += f' ({__git_sha__})'
