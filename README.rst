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

PyHDX is python project which can be used to derive Gibbs free energy and Protection Factors from HDX-MS data.
Currently the project is functional but in beta. Please refer to docs/installation.rst for installation instructions.

Installation 
===============
Below are the installation instructions from source:

1. Download or ``git clone`` the master branch of the PyHDX repository

2. Create a ``conda`` environment

    ``conda create --name <name> python=3.8``

3. Activate conda environment

    ``conda activate <name>``

4. Install the dependencies
    
    ``conda install -c conda-forge pyhdx --only-deps``

5. Building wheels for the project
    
    ``python setup.py sdist bdist_wheel``

6. Installing the wheels (should be generated in the dist folder)
    
    ``pip install dist/PyHDX-version.whl``

7. After successful installation, PyHDX web application can be launched from the command line using ``pyhdx`` command with below options,

    - To run PyHDX server using default settings on your local server, 
    
        ``pyhdx serve`` 
        
    - To run PyHDX server using the IP address and port number of your dask cluster, 
        
        ``pyhdx --cluster <ip> <port>``


Web Application
===============

A beta version of the web application is available for testing:
http://pyhdx.jhsmit.org/main

A test file can be downloaded from `here <https://raw.githubusercontent.com/Jhsmit/PyHDX/master/tests/test_data/ecSecB_apo.csv>`_. (right click, save as)


Two other web applications are available.
To upload fitting results from the main application and vizualize: 
http://pyhdx.jhsmit.org/single

To upload multiple fitting result datasets and compare and vizualize:
http://pyhdx.jhsmit.org/diff
