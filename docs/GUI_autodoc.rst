.. _apidocs:

Web Application
===============

This page contains auto-generated docs for PyHDX' web application.

There are three applications available:

 - **Main Application**
    Fitting of HDX-MS datasets, classification, visualization and exporting data.

 - **Single Classification**
    Reload exported data from the main application for classification, visualization and exporting data.

 - **Binary Comparison**
    Reload multiple exported datasets from the main application and calculate differences between pairs of datasets.
    The resulting differences can again be classified, visualized and exported.


The functionality in each app can be controlled by `Controllers` which can be found in the left sidebar. The functionality
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

