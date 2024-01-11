# Web application

This section will describe a typical workflow of using the main web interface application. Detailed information on each
parameter can be found in the web application [reference docs](web_app_autodoc.md). The web application consists of
a sidebar with controls and input, divided into sections, and a main view area with graphs and visualization. We will
go through the functionality of the web interface per section.

## Settings

In this section general settings for handling HDX data in the web application can be changed.

In the `Drop first`:material-form-textbox: field the number of N-terminal residues for each peptides can be chosen which 
should be ignored when calculating the maximum uptake for each peptide and are considered to fully 
exchange back. Prolines are ignored by default as they do not have exchangeable amide hydrogens.

The field `Weight Exponent`:material-form-textbox: controls how individual peptides are weighted when
calculating residue-level weighted averaged RFU values. The weight for each peptide is equal to
$\frac{1}{n^k}$ where the $n$ is equal to the number of exchanging residues, and $k$ the
user-configureable weight exponent. Therefore, shorter peptides contribute more to the
averaging result, as they contain higher resolution information. The default value for the exponent is 1,
and increasing the value increases the relative weight of shorter peptides.

## Peptide Input

Use `Input Mode`:material-form-dropdown:
to switch between input modes. When selecting '_Database_' HDX-MS datasets can be downloaded from the HDX-MS database (currently hosted on GitHub [here](https://github.com/Jhsmit/HDX-MS-datasets/)) and directly loaded into PyHDX. 

You can load your own Peptide d-uptake data into the web application either in '_Batch_' or '_Manual_' mode. 

When using '_Batch_', multiple measurements can be added quickly through an `.yaml` file specification.
An example of an `.yaml` file can be found in `tests/test-data` on GitHub ([here](https://github.com/Jhsmit/PyHDX/blob/master/tests/test_data/input/data_states.yaml)).

When using '_Single_', each measurement and experimental metadata must be input manually.
Use the `Browse`:material-button-cursor: button to select peptide data files to upload. These should be 'peptide master tables' which
is **long format** data where each entry should at least have the entries of:

 - `start` (inclusive residue number at which the detected peptide starts, first residue = 1)
 - `stop` (inclusive residue number at which the detected peptide stops)
 - `sequence` (sequence of the peptide in one letter amino acid codes)
 - `exposure` (time of exposure to deuterated solution)
 - `uptake` (amount/mass (g/mol) of deuterium taken up by the peptide)
 - `state` (identifier to which 'state' the peptide/protein is in (ie ligands, experimental conditions)

Currently, the only data format accepted is exported 'state data' from Waters DynamX, which is .csv format. Exposure
time units is assumed to be minutes. Other data format support can be added on request (eg HDExaminer).

Multiple files can be selected; one file can contain the experimental peptides while another file has the 
fully deuterated control sample peptides. When using a single file, the peptides should be marked using the 
'state' field.

Choose which method of back-exchange correction to use. Options are either to use a fully deuterated (FD) sample or
to set a fixed back-exchange percentage for all peptides. The latter method should only be used if no FD sample is
available. A percentage to set here can be obtained by running a back-exchange control once on your setup.

When selecting `FD Sample`:material-form-dropdown:, use the fields `FD File`:material-form-dropdown:, `FD State`:material-form-dropdown: 
and `FD Exposure`:material-form-dropdown: to choose which peptides should be used as FD control. 
Note that these peptides will be matched to the ones in the experiment and peptides without control will not be included.

Choose the experimental peptides using the `Exp File`:material-form-dropdown:, `Exp State`:material-form-dropdown:,
`Exp Exposure`:material-form-dropdown: and use `Experiment Exposures`:material-form-select: to select which 
exposure times to include.

Next, specify the percentage of deuterium in the labelling solution in the field `Deuterium percentage`:material-form-textbox:. This
percentage should be as high as possible; typically >=90%.

Use the fields `Temperature (K)`:material-form-textbox: and `pH read`:material-form-textbox: to specify the temperature 
and pH at which the D-labelling was done. The pH is the value as read from the pH meter without any correction. Temperature
units are Kelvin.

The next fields `N term`:material-form-textbox: and `C term`:material-form-textbox: specify the residue number indices 
of the N-terminal and C-terminal residues, respectively. For the
N-terminal this value is typically equal to 1, but if N-terminal affinity tags are used for purification this might be a
negative number. The value specified should match with the residue indices used in the input .csv file. The C-term value
tells the software at which index the C-terminal of the protein is, as it is possible that the protein extends beyond the
last residue included in any peptide and as the C-term exhibits different intrinsic rates of exchanges this needs to be
taken into account. A sequence for the full protein (in the N-term to C-term range as specified) can be added to provide
additional sequence information, but this is optional.

Finally, specify a name for the measurement, by default equal to the `Experiment State`:material-form-dropdown: value and 
press `Add Measurement` to add the measurement with the current specifications. Repeat the process to add additional 
measurements, either starting by adding a new file or changing the selection on the current file. The 
`Download HDX spec`:material-button-cursor: can be used to download a `.yaml` file with the full state 
specification, and this file can then in future sessions be used when using `Batch` input mode, by setting 
`Input mode`:material-form-dropdown: to '_Batch_'.

Finally, press the button `Load dataset`:material-button-cursor: to parse and load the full dataset.

Datasets currently cannot be removed, if you want to remove datasets, press the browser 'refresh' button to start over.

## Coverage (figure)

The '**Coverage**' figure in the main application area rectangles show corresponding to the peptides of a single
timepoint. Peptides are only included if they are in both all the timepoints as well as in the fully deuterated control
sample.

By hovering the mouse over the peptides in the graph, more information is shown about each peptide:

* `peptide_id`: Index of the peptide per timepoint starting at the first peptide at 0
* `start, end`: Inclusive, exclusive interval of residue numbers in this peptide (Taking N-terminal residues into account)
* `RFU`: Relative fraction uptake of the peptide
* `D(corrected)`: Absolute D-uptake, corrected by FD control
* `sequence`: FASTA sequence of the peptide. Non-exchanging N-terminal residues marked as 'x' and prolines in lower case.

## D-Uptake fit

This controller can be used to perform a fit of residue-level D-Uptake, independently for each timepoint in each measurement. The advantage of
these residue-level D-Uptake values compared to residue-level RFU values is that the former is calculated by weighted averaging, which results 
in loss of residue resolution, while the former is calculated by a least-squares fitting procedure, thus higher resulution results can be obtained
depending on the peptide overlap. D-Uptake fits can be advantageous over ΔG fits when HDX kinetics are not in the EX2 regime and approximations 
made in the ΔG fit procecure with respect to HDX kinetics are thus not applicable. 

D-Uptake fit take back exchange into account (by using the Fully Deuterated control sample) as well as the D-percentage of the labelling solution. Fully 
exchanged residues will therefore have a D-uptake value equal to the D fraction in solution (ie 0.9 for 90% deuterium).

As fitting of residue-level D-uptake also suffers from the issues of non-identifyability or underdetermined systems (more parameters than datapoints),
typically fits are repeated for N times with random initially guesses and a smoothing regularization term is applied along the primary structure. The number of 
fitting repeats can be set with `Repeats`:material-form-textbox:. The checkbox `Bounds`:material-checkbox-outline: controls whether or not bounds are applied to the fit (beteween 0 and 1), and 
should typically be checked. The value at `R1`:material-form-textbox: controls the degree of smoothing along primary structure. A good starting value is 1 (For the SecB test data),
and different values should be tried to find an optimal value and depends on the size of the protein and the number and sizes of available peptides. Finally, before 
starting a fit, choose a name for the fit with current settings at `Fit name`:material-form-textbox: and click `Do Fitting`:material-button-cursor: to start the fitting process. 

When the fit finishes, the found mean of the D-uptake values of all repeats will be shown as a scatterplot under the tab 'D-uptake'. The scatterplot has two sets of errorbars,
which show the 5, 25, 75 and 95 percentiles of the repeats; gray errorbars are 5-95 percentile and black errorbars are 25-75 percentiles. 

As with RFUs, differential D-uptake values between protein states can be calculated under `Differential HDX`:material-expand-all:. Both D-uptake and ΔD-uptake values can be directly 
visualized on the tertiary structure using `Protein Control`:material-expand-all: and the '**Protein**' viewer. 
When calculating and interpreting ΔD-uptake values, users should be aware that artefactual differences can arise as the fitting can converge to different solutions if in
specific regions not enough peptide overlap is available or if peptide coverage dramatically differs between protein states. Therefore, sanity checking of differential HDX
results with input data (peptide D-uptake / RFU values) is recommended.  


## Initial Guesses

As a first step in the fitting procedure, initial guesses for the exchange kinetics need to be derived. This can be done
through two options (`Fitting model`:material-form-dropdown:): 'Half-life' (fast but less accurate), or 'Association' (slower but more accurate).


Using the 'Association' procedure is recommended. This model fits two time constants to the weighted-averaged uptake kinetics of
each residue. At `Lower bound`:material-form-textbox: and `Upper bound`:material-form-textbox: the bounds of these rate constants can be specified
but in most cases the autosuggested bounds are sufficient. The bounds can be changed per dataset by using the `Dataset`:material-form-dropdown:
field or for all datasets at the same time by ticking the `Global bounds`:material-checkbox-outline: checkbox.
Rarely issues might arise when the initial guess rates are close to the specified bounds at which point the bounds should be
moved to contain a larger interval. This can be checked by comparing the fitted rates _k1_ and *k2* 
(`File Export`:material-expand-all: :material-arrow-right: `Target dataset`:material-form-dropdown: :material-arrow-right: `rates`).
Both rates and associated amplitudes are converted to a single rate value used for initial guesses.
To calculate guesses, select the model in the drop-down menu, assign a name to these initial guesses and the press
'Calculate Guesses'. The fitting is done in the background. When the fitting is done, the obtained rate is shown in the main area in the
tab 'Rates'. Note that these rates are merely a guesstimate of HDX rates and these rates should not be used for any
interpretation whatsoever but should only function to provide the global fit with initial guesses.

## ΔG Fit

After the initial guesses are calculated we can move on to the global fit of the data. Details of on the fitting procedure
can be found in the PyHDX publication in [Analytical Chemistry](https://doi.org/10.1021/acs.analchem.1c02155).

At `Initial guess`:material-form-dropdown:, select which dataset to use for initial guesses (typically '_Guess_1_'). Both previous fits (ΔG values)
as well as estimated HX rates can be used as initial guesses. The initial guesses can be applied as 'One-to-one', where each protein state
gets initial guesses derived from that state, or 'One-to-many', where one protein state is use as initial guesses for all states.
Users can switch between both modes using `Guess mode`:material-form-dropdown:.

At `Fit mode`:material-form-dropdown:, users can choose either '_Batch_' or '_Single_' fitting. If only one dataset is loaded, only '_Single_' is
available. If '_Single_' is selected, PyHDX will fit ΔG values for each dataset individually using the specified settings.
In '_Batch_' mode all data enters the fitting process at the same time. This allows for the use of a second regularizer
between datasets. Note that when using '_Batch_' mode, the relative magnitudes of the Mean Squared error losses and
regularizer might be different, such that '_Batch_' fitting with ``r2`` at zero is not identical to '_Single_' fits.

The fields `Stop loss`:material-form-textbox: and `Stop patience`:material-form-textbox: control the fitting termination. If the loss improvement
is less than `Stop loss`:material-form-textbox: for `Stop patience`:material-form-textbox: epochs (fit iterations), the fitting will terminate.
`Learning rate`:material-form-textbox: controls the step size per epoch. For typical a dataset with 62 peptides over 6 timepoints, the
learning rate should be 50-100. Smaller datasets require larger learning rates and vice versa.

`Momentum`:material-form-textbox: and `Nesterov`:material-checkbox-outline: are advanced settings passed to the Pytorch `SGD` optimizer.

The maximum number of epochs or fit iterations is set in the field `Epochs`:material-form-textbox:.

Finally, the fields `Regualizer 1`:material-form-textbox: and `Regualizer 2`:material-form-textbox: control the magnitude of the regualizers. Please refer
to our [publication](https://doi.org/10.1021/acs.analchem.1c02155) for more details. In short, `r1` acts along consecutive residues and affects as a 'smoothing'
along the primary structure. Higher values give a more smoothed result. This prevents overfitting or helps avoid problems
related to the 'non-identifiability' issue where in unresolved (no residue-level overlap) regions the correct kinetic components
can be found (ΔGs of residues given correct choice of timepoints) but it cannot confidently be assigned to residues as
resolution is lacking. The regualizer `r1` biases the fit result towards the residue assignment choice with the lowest
variation along the primary structure. Typical values range from 0.01 to 0.5, depending on size of the input data.

`r2` acts between samples, minimizing variability between them. This is used in differential HDX where users are interested
in ΔG differences (ΔΔG). When measuring HD exchange with differing experimental conditions, such as differences in peptides detected, timepoints
used or D-labelling temperature and pH, the datasets obtained will have different resolution, both 'spatially' (degree of
resolved residues) and 'temporally' (range/accuracy of ΔGs). This can lead to artefactual differences in the final ΔΔG result, as
features might be resolved in out dataset and not in the other, which will show up as ΔΔG.
The penalty from `r2` can be calculated either with respect to a selected reference state (Using `R2 reference`:material-form-dropdown:), or 
the average value between all states (set to `None`)

Specify a unique name at `Fit name`:material-form-textbox: and press `Do Fitting`:material-button-cursor: do start the fit. 
The '**Info log**' window in the bottom right corner displays information on when the fit started and finished. The fitting runs 
in the background and multiple jobs can be executed at the same time when processing multiple protein states with `Fit mode`:material-form-dropdown:
 set to 'Single'.
However, please take into account that these fits are computationally intensive and currently if multiple users submit too many jobs it might overwhelm our/your server.

The output ΔG values are shown in the '**ΔG**' graph.

See also the [fitting example notebook](examples/03_fitting.ipynb) for more details on fitting and the effect of regualizers.
