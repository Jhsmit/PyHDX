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


From sources
------------

1. Clone the PyHDX repository and cd into the project's root directory:
    .. code-block:: rst

        git clone https://github.com/Jhsmit/PyHDX.git
        cd PyHDX


2. Create a ``conda`` environment

    .. code-block:: rst

        conda create --name <name> python=3.8

3. Activate conda environment

    .. code-block:: rst

        conda activate <name>

4. Install PyTorch

If you would like a specific PyTorch version to use with PyHDX (ie CUDA/ROCm support), you should install this first.
Installation instructions are on the Pytorch_ website.

5. Install other dependencies

    .. code-block:: rst

        conda install -c conda-forge pyhdx=0.4.0b7 --only-deps``

This install dependencies only for base PyHDX. To install web application dependencies, 
run the `_requirements.py` file in the PyHDX root folder. This generates `reqs-<extra>.txt` files which lists
requirements.

    .. code-block:: rst
        python _requirements.py
        pip install -r reqs-base.txt -r req-web.txt -r req-web.txt

Or

    .. code-block:: rst
        python _requirements.py
        conda install --file reqs-base.txt --file req-web.txt --file req-web.txt


6. Install PyHDX in editable / development mode

    .. code-block:: rst

        conda develop .

    .. code-block:: rst

        pip install -e .

Dependencies
------------

The requirements for PyHDX are listed in setup.cfg

.. _Github repo: https://github.com/Jhsmit/pyhdx

.. _Pytorch: https://pytorch.org/
