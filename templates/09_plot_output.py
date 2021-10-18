#%%
from pyhdx.fileIO import load_fitresult
from pyhdx.plot import plot_fitresults
from pathlib import Path
import proplot as pplt
import matplotlib.pyplot as plt
import pandas as pd


#%%

# __file__ = Path().cwd() / 'templates'/ 'script.py'  # Uncomment for PyCharm scientific mode


cwd = Path(__file__).parent
output_dir = cwd / 'output' / 'figure'
fit_result = load_fitresult(cwd / 'output' / 'SecB_tetramer_dimer_batch')


plot_fitresults(cwd / 'output' / 'SecB_tetramer_dimer_batch')