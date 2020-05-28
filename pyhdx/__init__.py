# -*- coding: utf-8 -*-

"""Top-level package for PyHDX."""

__author__ = """Jochem Smit"""
__email__ = 'jhsmit@gmail.com'
__version__ = '0.1.0'

from .pyhdx import PeptideMasterTable, PeptideMeasurements, KineticsSeries, Coverage
from .fitting import KineticsFitting
from .fileIO import read_dynamx
from .output import Output
