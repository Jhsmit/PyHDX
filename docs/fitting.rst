=======
Fitting
=======

The main feature of pyHDX is the fitting of rate equations describing deuterium uptake to a kinetic series of measured
peptides each covering a section of residues with a corresponding amount of deuterium uptake per peptide-timepoint.

Overfitting
-----------

Overfitting occurs when more parameters are added to the model but the supplied data has insufficient independent datapoints
to be able to accurately and uniquely determine the value of these parameters. Typical signs of overfitting are large
variations along residues in the obtained rates, such as for residue 43 in Figure :numref:`overfitting`.

.. _overfitting:
.. figure:: figures/overfitting.png
    :scale: 25 %
    :figclass: align-center

    XX not really a great example of overfitting


To determine if overfitting has occurs, the number of fitting parameters should be varied while checking the effect of adding and removing fit parameters againts
goodness-of-fit parameters. This is a laborious and time consuming process and further streamlining and automating this
process is planned to be part of a future release.

In the current implementation, fitting accuracy and residue resolution is sacrificed in order to make sure overfitting is
unlikely. Block size is increased and the number of exchange rate time constants is limited to 2. The downside of this
approach is that the fits can be poor in the case of residues exchanging with more than two distinct rate constants per
block, or that features consisting of only several residues can be missed. Examples of how to customize the defintion of
fitting blocks can be found in the examples section.

Non-identifyability
-------------------

Consider a block of 5 amino acids which all exchange deuterium with very distinct exchange rates and a set of measurements
where the timepoints sufficiently cover these exchange rates. In this scenario, although its possible to extract all 5
kinetic rates by fitting the uptake curve, it is impossible to assign these kinetics rates to individual amino acids. This
is referred to as the non-identifyability issue (XX REF) and this can only be overcome by increasing the number of peptides
such that each amino acid occurs in a unique set of peptides.




