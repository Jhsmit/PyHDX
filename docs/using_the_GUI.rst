Using the web application
=========================


This section will describe a typical workflow of using the main web interface application. Detailed information on each
parameter can be found in the web application reference docs :ref:`web-application-autodoc`. The web application consists of
a sidebar with controls and input, divided into sections, and a main view area with graphs and visualization. We will
go through the functionality of the web interface per section.

Peptide Input
`````````````

Use the `Browse` button to select peptide data files to upload. These should be 'peptide master tables' which is **long format** data
where each entry should at least have the entries of:

 - start (inclusive residue number at which the detected peptide starts, first residue = 1)
 - stop (inclusive residue number at which the detected peptide stops)
 - sequence (sequence of the peptide in one letter amino acid codes)
 - exposure (time of exposure to deuterated solution)
 - uptake (amount/mass (g/mol) of deuterium taken up by the peptide)
 - state (identifier to which 'state' the peptide/protein is in (ie ligands, experimental conditions)

Currently the only data format accepted is exported 'state data' from Waters DynamX, which is .csv format. DynamX exposure
time units is assumed to be minutes. Other data format support can be added on request (eg HDExaminer).

Multiple files can be uploaded by using the 'Add File' button after which these files will be combined. Make sure there
are no overlaps/clashes between 'state' entries when combining multiple files.

In the 'Drop first' entry the number of N-terminal residues for each peptides can be chosen which should be ignored when
calculating the maximum uptake for each peptide as they are considered to fully exchange back. Prolines are ignored by
default as they do not have exchangeable amide hydrogens. Next, specify the percentage of deuterium in the labelling
solution. Then press 'Load Files'.

Next we need to specify how to correct for back-exchange. The recommended procedure is to use a fully deuterated (FD)
sample prepared by incubation in the deuterated solution for extended time periods or under denaturating conditions.
To do so, select the 'state' and 'exposure' fields of your FD sample under 'FD State' and 'FD Exposure'. Alternatively,
it is possible to set a fixed percentage of back-exchange by switching from experimental (Exp) to theoretical (Theory).

Now under 'Experiment State' and 'Experiment Exposures' choose which state and timepoints to process.

As the C-terminal residue in the protein experiences different intrinsic uptake kinetics, it is important to indicate which
residue index is the C-terminal residue. This should be the C-terminal index of the protein as used in the experiment, so
if C-terminal pufification tags are used this should be included.

Press 'Parse' to load the data and apply back-exchange corrections.

Coverage
````````
In the 'Coverage' figure in the main application area rectangles should show corresponding to the peptides of a single
timepoint. Peptides are only included if they are in both all the timepoints as well as in the fully deuterated control
sample. Under 'Coverage' in the sidebar the coverage graph can be controlled by specifying how many peptides to plot
vertically, which color map to use, which timepoint to show (using the slider) and which timepoint (Exposure) is
currently shown.

By hovering the mouse over the peptides in the graph, more information is shown about each peptide:

- Pos: current x (residue) position of the mouse
- Index: Index of the peptide per timepoint starting at the first peptide at 0
- Start: Inclusive index of the starting point of the peptide taking prolines and N-terminal residues into account. Original number from the input data is in brackets.
- End: Exclusive index of the end of the peptide, original number from the input data in brackets.
- Sequence: Sequence of the peptide with back-exchanging N-terminal residues marked as 'x' and prolines in lower case.
- Score: Percentage of deuterium uptake with respect to the maximum uptake as calculated from FD control, N-terminal residues and prolines.
- Uptake: Absolute amount of deuterium uptake as measured (Corrected uptake / number of exhangeing residues, max uptake as specified in the input data file).

Initial Guesses
```````````````

As a first step in the fitting procedure, initial guesses for the exchange kinetics need to be derived. This can be done
through two options: 'Half-life' (fast but less accurate), or 'Association' (slower but more accurate). Using the
'Association' procedure is recommended. This model fit two time constants the the weighted-averaged uptake kinetics of
each residue. A upper and lower bound of these rate constants can be specified but in most cases the autosuggested bounds
are sufficient.
Rarely issues might arise when the initial guess rates are close to the specified bounds at which point the bounds should be
moved. This can be checked by comparing the fitted rates *k1* and *k2* (File Export > export 'fit1') to the specified bounds.
Both rates are and associated amplitudes are converted to a single rate value used for initial guesses.
Select the model in the drop-down menu and the press 'Do fitting' to start fitting.
The fitting is done in the background a progress bar will show time until completion. When the fitting is done, the
obtained rate is shown in the main area in the tab 'Rates'. Currently the graph does not autoscale correctly so in order
to view the full data please use the navigation/zoom buttons on the right side of the graph.

Fitting
```````

After the initial guesses are calculated we can move on the the global fit of the data. Details of the fitting equation
can be found the PyHDX publication (currently `bioRxiv`_).

At 'Initial guess', select which dataset to use for initial guesses (typically 'fit1'). Then enter the temperature (in Kelvin)
of the labelling reaction and the pH of the labelling reaction (uncorrected pH read value).

The value of the regalizer control the degree of 'smoothing' which prevents overfitting. Typical values are 0.5 to 2, depending
on the input data, where lower values give more detail but should only be selected if the degree of peptide coverage and
overlap is high. For the other fitting hyperparameters, see the reference docs :ref:`web-application-autodoc`.

The output of the fit is ΔG, protection factor (PF), covariance (for ΔG) for each residue. All values can be exported in .txt
format and the ΔG and PF values are plotted in their respective graph windows.

Fit Results
```````````

The fit results panel controls the fit results graph where each peptide can be selected and measured (and corrected)
deuterium uptake values are plotted (with the fitted result, currently broken)

Classification
``````````````

The classification value can be used to calculate color assignments per residue from values of all available datasets.
Typically, ΔG values are used for classification. To do so, select 'global_fit' under 'Target' and 'deltaG' for 'quantity'.
This will calculate colors for the 'global_fit' dataset, if another column is subsequently used for coloring, for example
'pfact' (PF) or 'covariance' the colors are overwritten.

Two distince modes can be selected, 'Discrete', where all colors in a single defined category are the same, or 'Continuous',
where colors are interpolated linearly between defined nodes. This means that when three colors are chosen in the 'Discrete'
mode, two thresholds are defined to seperate the three classes, whereas in 'Continous' the number of thresholds is equal
to the number of colors.

The button 'Otsu' (only available in 'Discrete' coloring) automatically classifies values in the number of chosen categories
using Otsu's method (minimize variance within populations). With 'Linear' the thresholds are automatically equidistantly
spaced between the minimum and maximum value.

When the tickbox 'Log space' is selected this 'Linear' assignment is done in log space, as well as the color interpolation.
The thresholds as well as colors can be manually chosen. Note that the thresholds must always be decreasing in value from
Threshold 1.

File Export
```````````

The assigned colors per dataset as well as all datasets can be downloaded from the 'File Export' panel. Select the target
dataset to export and click the <name>_linear.txt button to export the raw data. For datasets which have an residue number
index column (r_number) have an additional pymol download button from which a .pml script can be downloaded. This script
can be ran from pymol to apply the colors to a 3D structure.

.. comment: check how the no coverage color is defined

Protein Viewer
``````````````

Assigned colors on a 3D structure can not only be exported to pymol but also directly visualized in the web application by
using the built in `NGL`_ protein viewer. A datasets should be selected which as previously assigned a color scheme in
**Classification**. Two structure input options are available, either a direct transfer from the RCSB PDB (choose Rcsb id
in the field below) or uploading a .pdb file.


.. _NGL: https://nglviewer.org
.. _bioRxiv: https://doi.org/10.1101/2020.09.30.320887