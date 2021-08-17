.. highlight:: shell

============
Installation
============


Stable release (v0.3.2)
-----------------------

Installation with `conda`:

.. code-block:: rst

   $ conda install -c conda-forge pyhdx

Installation with `pip`:

.. code-block:: rst

   $ pip install pyhdx==0.3.2


Beta release (v0.4.0b1)
-----------------------

To install base PyHDX:

.. code-block:: rst

   $ pip install pyhdx==0.4.0b1

To install with web application:

.. code-block:: rst

    $ pip install pyhdx==0.4.0b1[web]

Currently custom bokeh extensions are not packaged, therefore to run the web application Node.js is required:

.. code-block:: rst

    $ conda install nodejs

To install with pdf output:

.. code-block:: rst

    $ pip install pyhdx==0.4.0b1[pdf]

..
    From sources
    ------------

    1. Download or ``git clone`` the master branch of the PyHDX repository

    2. Create a ``conda`` environment

    .. code-block:: rst

        conda create --name <name> python=3.8

    3. Activate conda environment

    .. code-block:: rst

        conda activate <name>

    4. Install the dependencies

        ``conda install -c conda-forge pyhdx --only-deps``

    5. Building wheels for the project

        ``python setup.py sdist bdist_wheel``

    6. Installing the wheels (should be generated in the dist folder)

    ``pip install dist/PyHDX-version.whl``


Running the web server
----------------------

PyHDX web application can be launched from the command line using ``pyhdx`` command with below options,

To run PyHDX server using default settings on your local server:

.. code-block:: rst

    $ pyhdx serve

To run PyHDX server using the IP address and port number of your dask cluster:

.. code-block:: rst

    $ pyhdx --scheduler_address <ip>:<port>

If no dask cluster is found at the specified address, a LocalCluster will be started (on localhost) using the
specified port number.

To start a dask cluster separately, open another terminal tab and run:

.. code-block:: rst

    python local_cluster.py

This will start a Dask cluster on the scheduler address as specified in the PyHDX config.
(user dir / .pyhdx folder)


Dependencies
------------

The requirements for PyHDX are listed in setup.cfg

.. _Github repo: https://github.com/Jhsmit/pyhdx
