# Adapted from: https://github.com/maresb/hatch-vcs-footgun-example
# Define the variable '__version__':
try:
    # If we are in an editable install, the _versioneer file exist and we can use it to find the version
    from pyhdx._versioneer import get_versions

    # This will fail with LookupError if the package is not installed in
    # editable mode or if Git is not installed.
    __version__ = get_versions()["version"]
except ImportError:
    # If the project build with hatch, there should be a _version.py file
    try:
        from pyhdx._version import __version__  # noqa: F401 # type: ignore
    except ModuleNotFoundError:
        # The user is probably trying to run this without having installed
        # the package, so complain.
        raise RuntimeError("PyHDX is not correctly installed. Please install it with pip.")
