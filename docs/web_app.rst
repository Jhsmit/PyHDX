PyHDX web application
=====================

This section will describe a typical workflow of using the main web interface application. Detailed information on each
parameter can be found in the web application reference docs :ref:`web-application-autodoc`. The web application consists of
a sidebar with controls and input, divided into sections, and a main view area with graphs and visualization. We will
go through the functionality of the web interface per section.

Peptide Input
`````````````

Use the :guilabel:`Browse` button to select peptide data files to upload. These should be 'peptide master tables' which
is **long format** data
where each entry should at least have the entries of:

 - start (inclusive residue number at which the detected peptide starts, first residue = 1)
 - stop (inclusive residue number at which the detected peptide stops)
 - sequence (sequence of the peptide in one letter amino acid codes)
 - exposure (time of exposure to deuterated solution)
 - uptake (amount/mass (g/mol) of deuterium taken up by the peptide)
 - state (identifier to which 'state' the peptide/protein is in (ie ligands, experimental conditions)

Currently the only data format accepted is exported 'state data' from Waters DynamX, which is .csv format. Exposure
time units is assumed to be minutes. Other data format support can be added on request (eg HDExaminer).

Multiple files can be selected after which these files will be combined. Make sure there are no overlaps/clashes
between 'state' entries when combining multiple files.

Choose which method of back-exchange correction to use. Options are either to use using a fully deuterated sample or
to set a fixed back-exchange percentage for all peptides. The latter method should only be used if no FD sample is
available. A percentage to set here can be obtained by running a back-exchange control once on your setup.

When selecting :guilabel:`FD Sample`, use the fields :guilabel:`FD State` and :guilabel:`FD Exposure` to choose which
peptides from the input should be used as FD control. Note that these peptides will be matched to the ones in the
experiment and peptides without control will not be included.

Use the fields :guilabel:`Experiment State` to choose the 'state' of your experiment. In 'Experiment Exposures' you can select
which exposure times to add include the dataset.

In the :guilabel:`Drop first`'` entry the number of N-terminal residues for each peptides can be chosen which should be ignored when
calculating the maximum uptake for each peptide as they are considered to fully exchange back. Prolines are ignored by
default as they do not have exchangeable amide hydrogens.

Next, specify the percentage of deuterium in the labelling solution in the field :guilabel:`Deuterium percentage`. This
percentage should be as high as possible, typically >90%.

Use the fields :guilabel:`Temperature (K)` and :guilabel:`pH read` to specify the temperature and pH at which the D-labelling
was done. The pH is the value as read from the pH meter without any correction.

The next fields :guilabel:`N term` and :guilabel:`C term` specify the residue number indices of the N-terminal and C-terminal residues, respectively. For the
N-terminal this value is typically equal to 1, but if N-terminal affinity tags are used for purification this might be a
negative number. The value specified should match with the residue indices used in in the input .csv file. The C-term value
tells the software at which index the C-terminal of the protein is, as it is possible that the protein extends beyond the
last residue included in any peptide and as the C-term exhibits different intrinsic rates of exchanges this needs to be
taken into account. A sequence for the full protein (in the N-term to C-term range as specified) can be added to provide
additional sequence information, but this is optional.

Finally, specify a name of the dataset, by default equal to the 'state' value and press 'Add dataset' to add the dataset.
Datasets currently cannot be removed, if you want to remove datasets, press the browser 'refresh' button to start over.

Coverage
````````

In the 'Coverage' figure in the main application area rectangles should show corresponding to the peptides of a single
timepoint. Peptides are only included if they are in both all the timepoints as well as in the fully deuterated control
sample. Under 'Coverage' in the sidebar, part of the coverage graph can be controlled, by specifying which color map to use.
Controls for choosing which 'state' to display and which exposure time can be found under Graph Cont

..
    #how many peptides to plot
    vertically, which color map to use, which timepoint to show (using the slider) and which timepoint (Exposure) is
    currently shown.

By hovering the mouse over the peptides in the graph, more information is shown about each peptide:

* index: Index of the peptide per timepoint starting at the first peptide at 0
* rfu: Relative fraction uptake of the peptide (Absolute corrected D-uptake in brackets)
* sequence: FASTA sequence of the peptide. Non-exchanging N-terminal reisues marked as 'x' and prolines in lower case.

..
    - Start: Inclusive index of the starting point of the peptide taking prolines and N-terminal residues into account. Original number from the input data is in brackets.
    - End: Exclusive index of the end of the peptide, original number from the input data in brackets.

Initial Guesses
```````````````

As a first step in the fitting procedure, initial guesses for the exchange kinetics need to be derived. This can be done
through two options (:guilabel:`Fitting model`): 'Half-life' (fast but less accurate), or 'Association' (slower but more accurate). Using the
'Association' procedure is recommended. This model fits two time constants the the weighted-averaged uptake kinetics of
each residue. At :guilabel:`Lower bound` and :guilabel:`Upper bound` the bounds of these rate constants can be specified
but in most cases the autosuggested bounds are sufficient. The bounds can be changed per dataset by using the :guilabel:`Dataset`
field or for all datasets at the same time by ticking the :guilabel:`Global bounds` checkbox.
Rarely issues might arise when the initial guess rates are close to the specified bounds at which point the bounds should be
moved to contain a larger interval. This can be checked by comparing the fitted rates *k1* and *k2* (:menuselection:`File Export --> Target dataset --> rates`)
Both rates and associated amplitudes are converted to a single rate value used for initial guesses.
To calcualte guesses, select the model in the drop-down menu, assign a name to these initial guesses and the press
'Calculate Guesses'. The fitting is done in the background. When the fitting is done, the obtained rate is shown in the main area in the
tab 'Rates'. Note that these rates are merely an guesstimate of HDX rates and these rates should not be used for any
interpretation whatsoever but should only function to provide the global fit with initial guesses.

Fitting
```````

After the initial guesses are calculated we can move on the the global fit of the data. Details of the fitting equation
can be found the PyHDX publication (currently `bioRxiv`_).

At 'Initial guess', select which dataset to use for initial guesses (typically 'Guess_1').
At 'Fit mode', users can choose either 'Batch' or 'Single' fitting. If only one datasets is loaded, only 'Single' is
available. If 'Single' is selected, PyHDX will fit ΔG values for each datasets individually using the specified settings.
In 'Batch' mode all data enters the fitting process at the same time. This allows for the use of a second regularizer
between datasets. Note that when using 'Batch' mode, the relative magnitudes of the Mean Squared error losses and
regularizer might be different, such that 'Batch' fitting with ``r2`` at zero is not identical to 'Single' fits.

The fields :guilabel:`Stop loss` and :guilabel:`Stop patience` control the fitting termination. If the loss improvement
is less than `Stop loss` for `Stop patience` epochs (fit iterations), the fitting will terminate.
:guilabel:`Learning rate` controls the step size per epoch. For typical a dataset with 62 peptides over 6 timepoints, the
learning rate should be 50-100. Smaller datasets require larger learning rates and vice versa.

:guilabel:`Momentum` and :guilabel:`Nesterov` are advanced settings for the Pytorch ``SGD`` optimizer.

The maximum number of epochs or fit iterations is set in the field :guilabel:`Epochs`.

Finally, the fields :guilabel:`Regualizer 1` and :guilabel:`Regualizer 2` control the magnitude of the regualizers. Please refer
to our `bioRxiv`_ manuscript for more details. In short, ``r1`` acts along consecutive residues and affects as a 'smoothing'
along the primary structure. Higher values give a more smoothed result. This prevents overfitting or helps avoid problems
in the 'non-identifiability' issue where in unresolved (no residue-level overlap) regions the correct kinetic components
can be found (ΔGs of residues given correct choice of timepoints) but it cannot confidently be assigned to residues as
resolution is lacking. The regualizer `r1` biases the fit result towards the residue assignment choice with the lowest
variation along the primary structure. Typical values range from 0.01 to 0.5, depending on size of the input data.

`r2` acts between samples, minimizing variability between them. This is used in differential HDX where users are interested
in ΔG differences (ΔΔG). When measuring HD exchange with differing experimental conditions, such as differences in peptides detected, timepoints
used or D-labelling temperature and pH, the datasets obtained will have different resolution, both 'spatially' (degree of
resolved residues) and 'temporally' (range/accuracy of ΔGs). This can lead to artefactual differences in the final ΔΔG result, as
features might be resolved in out dataset and not in the other, which will show up as ΔΔG.

Specify a unique name at :guilabel:`Fit name` and press :guilabel:`Do Fitting` do start the fit. The :guilabel:`Info log`
in the bottom right corner displays information on when the fit started and finished. The fitting runs in the background
and multiple jobs can be executed at the same time. However, please take into account that these fits are computationally
intensive and currently if multiple users submit too many jobs it might overwhelm our/your server.

The output ΔG values are shown in the 'Gibbs' graph (bottom left).

See also the :doc:`Fitting example <../examples/03_fitting>` section for more details on fitting and the effect of regualizers.

Graph Control
`````````````

This section is used to control which dataset is currently show in the graphs. Use the selector :guilabel:`Fit id` to
switch between fit results. The selector :guilabel:`State name` is used to switch between experimental states.
In the current version (v0.4.0bx) the switching is a bit slow so please wait a few seconds before switching again.

Under :guilabel:`Coverage`, use the :guilabel:`Exposure` selector to choose which timepoint is shown in the the coverage
graph. At :guilabel:`Peptide` and :guilabel:`Peptide index` the peptide to show at the Peptide (uptake curve) graph can
be chosen.

We can use these control to inspect the quality of the fit obtained. First, at :guilabel:`Losses` (bottom right) the progress
of the fit can be inspected. This should show a rapid decrease of the 'mse' loss, followed by a mostly flat plateau. If this
is not the case, extend the number of epochs (:guilabel:`epochs` or :guilabel:`stop_loss` and :guilabel:`Stop patience`)
or increase :guilabel:`Learning rate`.

The graph 'coverage mse' shows the total mean squared error of all timepoints per peptide. The color scale adjust automatically
so red colors do not necessarily reflect a poor fit, but highligh the worst fitted peptides in your dataset. Hover over
the peptide with the mouse to find the index of the peptide and select the peptide with :guilabel:`Peptide index`.

Classification
``````````````

The classification value can be used to create color datasets based on results from the global fit. At :guilabel:`fit_ID`,
choose which of the fit runs to use. Use :guilabel:`state_name` to choose which experimental states to apply the
color map to, use '`*`' to select all states. Finally use :guilabel:`quantity` to select which output column to use (typically
deltaG)

A preview of the color map will be applied to the ΔG values in the 'Gibbs' graph.

:guilabel:`Mode` can be used to select between the available color modes; `Discrete`, `Continuous` and `Color map`. `Discrete`
splits the ΔG values in `n` categories, which are all assigned the same color. When using `Continuous`, `n` color 'nodes' can be
defined, where color values are interpolated between these nodes. `Color map` allows users to choose a colormap from either
``matplotlib`` or ``colorcet``.

The number of categories can be set with :guilabel:`Number of colours`.
When using `Discrete` coloring, the thresholds of the categories can be automatically determined by pressing the :guilabel:`Otsu`
button (using Otsu's method). Use the button :guilabel:`Linear` to distribute threshold values automatically with equal
distances between them, and the extrema at the largest/smallest ΔG values.

Toggle :guilabel:`Log space` to apply the color map in log space (typically used for colouring rates/protection factors).
Assign an unique name using :guilabel:`Color set name` and press :guilabel:`Add colorset` to create the color dataset.

A color for residues which are covered by peptides can be chosen at :guilabel:`No coverage`.
The colors for the color groups or nodes can be chosen at the bottom of the controllers, as well as the exact position
of the thresholds. These values must be input such that they are always in decreasing order.

Protein Viewer
``````````````

Colorsets can be directly visualized on a protein structure using the built in `NGL`_ protein viewer. Use :guilabel:`Input mode`
to choose whether to download a structure from the RCSB database or to use a local file. When selecting `RCSB`, enter
the entry ID under :guilabel:`Rcsb id` and press :guilabel:`Load structure` to download the pdb file and display it. When
selecting `PDB File` use the :guilabel:`Choose File` button to upload a .pdb file from your computer, and press :guilabel:`Load structure`

Select wich color map to display on the structure with :guilabel:`color ID`, and select the experimental state with
:guilabel:`State name`

File Export
```````````

All tables which underlie the graphs in the PyHDX web application can be downloaded directly. Choose the the desired dataset
('global_fit' for ΔG values) with :guilabel:`Target dataset`. The data can be exported in machine-readable .csv files or
human-readable .txt (pprint) file by setting :guilabel:`Export format`. Make sure to download at least the .csv file for
further.

When selecting a color dataset, the data can not only be download as a .csv file but also as (a zip file of) .pml files
which contain pymol scripts to directly apply the color map to a structure in pymol.


.. _NGL: https://nglviewer.org
.. _bioRxiv: https://doi.org/10.1101/2020.09.30.320887
