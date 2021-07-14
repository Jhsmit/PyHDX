"""The setup script."""

from setuptools import setup
import versioneer

test_requirements = ['pytest', ]

setup(
    version=versioneer.get_version(),
)
