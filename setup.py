"""The setup script."""

from setuptools import setup
import versioneer


test_requirements = ['pytest', ]

#setup_requires=['pbr>1.9', 'setuptools>17.1'],
setup(
    version=versioneer.get_version(),
    packages=["pyhdx"],
    entry_points={
         "console_scripts": [
             "pyhdx = pyhdx.cli:main",

         ],
     }
)
