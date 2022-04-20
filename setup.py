"""The setup script."""

from setuptools import setup, find_packages
import versioneer

setup(
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    include_package_data=True,
)
