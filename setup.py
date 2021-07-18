"""The setup script."""

from setuptools import setup, find_packages
import versioneer

setup(
    version=versioneer.get_version(),
    #version='v0.4.0-rc3',
    packages=find_packages(),
    include_package_data=True
)
