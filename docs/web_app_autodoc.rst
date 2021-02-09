.. _web-application-autodoc:

Web Application Reference
=========================

This page contains auto-generated docs for PyHDX' web application.

There are five applications available:

- **Main Application** (main)
    Fitting of HDX-MS datasets (PFs and ΔG), classification, visualization and exporting data.

- **Single Classification** (single)
    Reload exported data from the main application for classification, visualization and exporting data.

- **Binary Comparison** (diff)
    Reload multiple exported datasets from the main application and calculate differences between pairs of datasets.
    The resulting differences can again be classified, visualized and exported.

- **Full deuteration** (full_deuteration)
    Load a fully deuterated dataset and check the amount of deuteration as a percentage of the theoretical maximum amount of deuterium uptake.

- **Folding** (folding)
    Highly experimental app for the analysis of folding data.

The functionality in each app can be controlled by `Controllers` which can be found in the left sidebar. The control parameters
of every controller per app is listed in the sections below.


Main Application
^^^^^^^^^^^^^^^^

.. autoclass:: pyhdx.panel.controllers.PeptideFileInputControl
.. autoclass:: pyhdx.panel.controllers.CoverageControl
.. autoclass:: pyhdx.panel.controllers.InitialGuessControl
.. autoclass:: pyhdx.panel.controllers.FitControl
.. autoclass:: pyhdx.panel.controllers.ClassificationControl
.. autoclass:: pyhdx.panel.controllers.FileExportControl
.. autoclass:: pyhdx.panel.controllers.ProteinViewControl
.. autoclass:: pyhdx.panel.controllers.OptionsControl

Single Classification
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: pyhdx.panel.controllers.MappingFileInputControl
.. autoclass:: pyhdx.panel.controllers.SingleControl
.. autoclass:: pyhdx.panel.controllers.ClassificationControl
.. autoclass:: pyhdx.panel.controllers.ProteinViewControl
.. autoclass:: pyhdx.panel.controllers.DifferenceFileExportControl
.. autoclass:: pyhdx.panel.controllers.OptionsControl

Binary Comparison
^^^^^^^^^^^^^^^^^

.. autoclass:: pyhdx.panel.controllers.MappingFileInputControl
.. autoclass:: pyhdx.panel.controllers.DifferenceControl
.. autoclass:: pyhdx.panel.controllers.ClassificationControl
.. autoclass:: pyhdx.panel.controllers.ProteinViewControl
.. autoclass:: pyhdx.panel.controllers.DifferenceFileExportControl
.. autoclass:: pyhdx.panel.controllers.OptionsControl

Full Deuteration
^^^^^^^^^^^^^^^^

.. autoclass:: pyhdx.panel.controllers.FDPeptideFileInputControl
.. autoclass:: pyhdx.panel.controllers.FDCoverageControl
.. autoclass:: pyhdx.panel.controllers.OptionsControl


Folding
^^^^^^^

.. autoclass:: pyhdx.panel.controllers.PeptideFoldingFileInputControl
.. autoclass:: pyhdx.panel.controllers.CoverageControl
.. autoclass:: pyhdx.panel.controllers.FoldingFitting
.. autoclass:: pyhdx.panel.controllers.FitResultControl
.. autoclass:: pyhdx.panel.controllers.ClassificationControl
.. autoclass:: pyhdx.panel.controllers.ProteinViewControl
.. autoclass:: pyhdx.panel.controllers.FileExportControl
.. autoclass:: pyhdx.panel.controllers.OptionsControl