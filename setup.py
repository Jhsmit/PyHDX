#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0', ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

#setup_requires=['pbr>1.9', 'setuptools>17.1'],
setup(pbr=True)

# setup(
#     author="Jochem Smit",
#     author_email='jhsmit@gmail.com',
#     classifiers=[
#         'Development Status :: 2 - Pre-Alpha',
#         'Intended Audience :: Developers',
#         'Natural Language :: English',
#         'Programming Language :: Python :: 3.6',
#         'Programming Language :: Python :: 3.7',
#     ],
#     description="PyHDX",
#     entry_points={
#         'console_scripts': [
#             'pyhdx=pyhdx.cli:main',
#         ],
#     },
#     install_requires=requirements,
#     long_description=readme + '\n\n' + history,
#     include_package_data=True,
#     keywords='pyhdx',
#     name='pyhdx',
#     packages=find_packages(include=['pyhdx']),
#     setup_requires=setup_requirements,
#     test_suite='tests',
#     tests_require=test_requirements,
#     url='https://github.com/Jhsmit/pyhdx',
#     version='0.1.0',
#     zip_safe=False,
# )
