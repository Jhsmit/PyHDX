.. highlight:: shell

============
Installation
============


Stable release
--------------

Installation with `conda`:

.. code-block::
    $ conda install -c conda-forge pyhdx

Installation with `pip`:

.. code-block::
    $ pip install pyhdx



From sources
------------

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


Dependencies
------------

The requirements for PyHDX are listed in requirements.txt

.. _Github repo: https://github.com/Jhsmit/pyhdx
