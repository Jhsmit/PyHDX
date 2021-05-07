=====
PyHDX
=====

|biorxiv| |test_action| |zenodo| |license| |docs| |coverage|

.. |zenodo| image:: https://zenodo.org/badge/206772076.svg
   :target: https://zenodo.org/badge/latestdoi/206772076

.. |biorxiv| image:: https://img.shields.io/badge/bioRxiv-v1-%23be2635
   :target: https://www.biorxiv.org/content/10.1101/2020.09.30.320887v1
   
.. |license| image:: https://img.shields.io/badge/License-MIT-yellow.svg
    :target: https://opensource.org/licenses/MIT

.. |test_action| image:: https://github.com/Jhsmit/PyHDX/workflows/pytest/badge.svg
    :target: https://github.com/Jhsmit/PyHDX/actions?query=workflow%3Apytest
    
.. |docs| image:: https://readthedocs.org/projects/pyhdx/badge/?version=latest
    :target: https://pyhdx.readthedocs.io/en/latest/?badge=latest

.. |coverage| image:: https://codecov.io/gh/Jhsmit/PyHDX/branch/master/graph/badge.svg?token=PUQAEMAUHH
      :target: https://codecov.io/gh/Jhsmit/PyHDX
    

`Documentation <https://pyhdx.readthedocs.io>`_

.. raw:: html

    <img src="images/PyHDX_rgb.png" width="1000" />

PyHDX is python project which can be used to derive Gibbs free energy from HDX-MS data.
Currently the project is functional but in beta, future versions will likely subject to changes in API and analysis.

Installation 
============

Installation with `conda`:

.. code-block::

    $ conda install -c conda-forge pyhdx

Installation with `pip`:

.. code-block::

    $ pip install pyhdx

Development install instructures in docs/installation.rst

Run PyHDX
=========

.. code-block::

    $ pyhdx serve
    
Please refer to the `docs <https://pyhdx.readthedocs.io>`_ for more details on how to run PyHDX.


Web Application
===============

The PyHDX web application is currently hosted at:
http://pyhdx.jhsmit.org/main

A test file can be downloaded from `here <https://raw.githubusercontent.com/Jhsmit/PyHDX/master/tests/test_data/ecSecB_apo.csv>`_. (right click, save as)


Two other web applications are available.
To upload fitting results from the main application and vizualize: 
http://pyhdx.jhsmit.org/single

To upload multiple fitting result datasets and compare and vizualize:
http://pyhdx.jhsmit.org/diff
