============
Installation
============

Currently the recommended version to use is the latest beta release (v0.4.0bx)

Stable release (v0.3.2)
-----------------------

Installation with `conda`:

.. code-block:: rst

   $ conda install -c conda-forge pyhdx

Installation with `pip`:

.. code-block:: rst

   $ pip install pyhdx==0.3.2


Beta release (v0.4.0b8)
-----------------------

To install base PyHDX:

.. code-block:: rst

   $ pip install pyhdx==0.4.0b8

To install with web application:

.. code-block:: rst

    $ pip install pyhdx==0.4.0b8[web]

To install with pdf output:

.. code-block:: rst

    $ pip install pyhdx==0.4.0b8[pdf]


Install from source
-------------------

Create a new conda environment:

.. code-block::

    $ conda create --name py38_pyhdx python=3.8
    # conda activate py38_pyhdx

Clone the github repository:

.. code-block:: rst

    $ git clone https://github.com/Jhsmit/PyHDX
    $ cd PyHDX

Generate conda requirements files from `setup.cfg`:

.. code-block:: rst

    $ python _requirements.py

Install the base dependencies and optional extras. For example, to install PyHDX with web app:

.. code-block:: rst

    $ conda install --file _req-base.txt --file _req-web.txt

To run the web application:

.. code-block::

    $ python pyhdx/web/serve.py

This runs the pyhx web application without a Dask cluster to submit jobs to, so
submitting a fitting job will give an error.

To start a dask cluster separately, open another terminal tab and run:

.. code-block:: rst

    python local_cluster.py



Running the web server
----------------------

PyHDX web application can be launched from the command line using ``pyhdx`` command with below options,

To run PyHDX server using default settings on your local server:

.. code-block:: rst

    $ pyhdx serve

To run PyHDX server using the IP address and port number of your dask cluster:

.. code-block:: rst

    $ pyhdx serve --scheduler_address <ip>:<port>

If no dask cluster is found at the specified address, a LocalCluster will be started (on localhost) using the
specified port number.

To start a dask cluster separately, open another terminal tab and run:

.. code-block:: rst

    python local_cluster.py

This will start a Dask cluster on the scheduler address as specified in the PyHDX config.
(user dir / .pyhdx folder)


Install from source
-------------------

Create a new conda environment:

.. code-block::

    $ conda create --name py38_pyhdx python=3.8
    # conda activate py38_pyhdx

Clone the github repository:

.. code-block:: rst

    $ git clone https://github.com/Jhsmit/PyHDX
    $ cd PyHDX

Generate conda requirements files `from setup.cfg`:

.. code-block:: rst

    $ python _requirements.py


If you would like a specific PyTorch version to use with PyHDX (ie CUDA/ROCm support), you should install this first.
Installation instructions are on the Pytorch_ website.

Then, install the other base dependencies and optional extras. For example, to install PyHDX with web app:

.. code-block:: rst

    $ conda install --file _req-base.txt --file _req-web.txt

Optionally, install PyHDX in develop/editable mode

.. code-block:: rst

    conda develop .

To run the web application:

.. code-block::

    $ python pyhdx/web/serve.py

This runs the pyhx web application without a Dask cluster to submit jobs to, so
submitting a fitting job will give an error.

To start a dask cluster separately, open another terminal tab and run:

.. code-block:: rst

    python local_cluster.py


Dependencies
------------

The requirements for PyHDX and its extras are listed in setup.cfg

.. _Github repo: https://github.com/Jhsmit/pyhdx

.. _Pytorch: https://pytorch.org/
