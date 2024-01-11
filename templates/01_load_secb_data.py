"""Load a HDX-MS dataset (peptide list with D-uptake per peptide, csv format)"""
# %%
from pathlib import Path
import numpy as np
from pyhdx import read_dynamx, HDXMeasurement
from pyhdx.process import apply_control, correct_d_uptake
from pyhdx.datasets import filter_peptides

# %%

current_dir = Path(__file__).parent
output_dir = current_dir / "output"
output_dir.mkdir(exist_ok=True)
np.random.seed(43)

# %%

fpath = current_dir.parent / "tests" / "test_data" / "input" / "ecSecB_apo.csv"

# Load the full .csv file
df = read_dynamx(fpath)

fd = {"state": "Full deuteration control", "exposure": {"value": 0.167, "unit": "min"}}

# Filter out peptides for the full deuteration control
fd_df = filter_peptides(df, **fd)

# filter out peptides for the experiment with state "SecB WT apo"
peptides = filter_peptides(df, state="SecB WT apo")

# Apply FD control, returns only peptides in both FD control and experiment
peptides_control = apply_control(peptides, fd_df)

# Correct for N-terminal back exchanging residues and deuterium percentage of the buffer
peptides_corrected = correct_d_uptake(peptides_control, drop_first=1, d_percentage=90.0)

sequence = "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA"
temperature, pH = 273.15 + 30, 8.0

# %%

# Create HDX Measurement object with addtional experimental metadata (sequence, pH, temperature)
hdxm = HDXMeasurement(peptides_corrected, sequence=sequence, pH=pH, temperature=temperature)

#%%
data = hdxm.data
data['id'] = data['start'].apply(str) + '_' + data['end'].apply(str)
peptides = data['id'].unique()
peptides
#%%

peptide = peptides[15]
peptide = '118_125'

df_f = data[data['id'] == peptide]
row_0 = df_f.iloc[0]
start, end = row_0['start'], row_0['end']
start, end

df_f

#%%
from hdxrate import k_int_from_sequence
seq = df_f['_sequence'].iloc[0]

k_int_from_sequence(seq, temperature, pH, d_percentage=90)

exposure = np.array(df_f['exposure'])[1:]
upt = np.array(df_f['uptake_corrected'])[1:]

if exposure[0] != 0:
    exposure = np.insert(exposure, 0, 0)
    upt = np.insert(upt, 0, 0)

exposure, upt

#%%
import proplot as pplt
fig, ax = pplt.subplots()
ax.scatter(exposure, upt)
ax.format(xscale='log')
pplt.show()

#%%
peptide_info = hdxm.coverage.protein.loc[start:end]
k_int = peptide_info['k_int'][peptide_info['exchanges']]
k_arr = np.array(k_int)

k_eff = np.exp(np.log(k_arr).mean())

k_int = np.array(peptide_info['k_int'])
k_int

# TODO i need the inverse of this function to find the time difference; then take the log of it
d0 = 1 - np.exp(np.divide(-k_int[:, np.newaxis]*exposure[np.newaxis, :], 2))
d0

#%%
np.log10(10)

#%%

np.log(2.72)

#%%

hdxm.metadata

#%%
 
data_wide = data.pivot(index='id', columns='exposure', values='uptake_corrected')

max_d = data[data['id'] == peptide]['ex_residues'].iloc[0]

import proplot as pplt

t_eval = np.linspace(0, 6000*10, 100, endpoint=True)
k = 1e-3
d_eval = max_d * (1 - np.exp(-k*t_eval))
#%%

from scipy.interpolate import PchipInterpolator

s = data_wide.loc[peptide]
x = np.array(s.index)
y = np.array(s)

x = np.append(x, x[-1]*5)
y = np.append(y, y[-1])

from scipy.interpolate import make_interp_spline

bspl = make_interp_spline(x, y, k=2)


#%%
# parameterized spline
p = np.stack((x, y))

u_unif = x

dp = p[:, 1:] - p[:, :-1]      # 2-vector distances between points
l = (dp**2).sum(axis=0)        # squares of lengths of 2-vectors between points
u_cord = np.sqrt(l).cumsum()   # cumulative sums of 2-norms
u_cord = np.r_[0, u_cord]      # the first point is parameterized at zero

u_c = np.r_[0, np.cumsum((dp**2).sum(axis=0)**0.25)]

from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 3, figsize=(8, 3))
parametrizations = ['uniform', 'cord length', 'centripetal']
for j, u in enumerate([u_unif, u_cord, u_c]):
   spl = make_interp_spline(u, p, axis=1)    # note p is a 2D array
   uu = np.linspace(u[0], u[-1], 51)
   xx, yy = spl(uu)
   ax[j].plot(xx, yy, '--')
   ax[j].plot(p[0, :], p[1, :], 'o')
   ax[j].set_title(parametrizations[j])

#plt.show()

xx, uu
#%%

spl_p = make_interp_spline(u_c, p, axis=1)  
uu = np.linspace(u_c[0], u_c[-1], 51)
xx, yy = spl_p(uu)

spl = PchipInterpolator(x, y)

fig, ax = pplt.subplots()
ax.scatter(data_wide.loc[peptide])
ax.plot(t_eval, d_eval, color='k', ls='--')
ax.plot(t_eval, spl(t_eval), color='g', ls='--')
#ax.plot(t_eval, bspl(t_eval), color='b', ls='--')
ax.plot(xx, yy, color='c', ls='--')
# ax.axhline(max_d, color='r', ls='--')
# ax.format(xscale='log')
ax.format(xlim=(0, 100), ylim=(0, 20))

#%%
exposures = np.array(data_wide.columns)
exposures

#%%

- np.log(1/2), np.log(2)

#%%
hdxm.data_wide.iloc[0].xs('d_uptake')

#%%

