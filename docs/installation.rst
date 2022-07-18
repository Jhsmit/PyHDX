============
Installation
============

Stable release
--------------

Installation of the latest stable release version.

With `conda`:

.. code-block:: rst

   $ conda install -c conda-forge pyhdx

With `pip`:

.. code-block:: rst

   $ pip install pyhdx

To install with web application:

.. code-block:: rst

    $ pip install pyhdx[web]

To install with pdf output:

.. code-block:: rst

    $ pip install pyhdx[pdf]


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

Dependencies can then be installed either from pinned versions or by using pip/conda to
resolve dependencies.

Pinned Dependencies
*******************

You can use one of the files in 'requirements/pinned' to install a pretested set of pinned
dependencies.

With `pip`:

.. code-block:: rst

    $ pip install -r requirements/pinned/py38_windows_pip.txt

Or use 'py38_linux_pip.txt' (These files should be the same)

With `conda`:

.. code-block:: rst

    $ conda env create -f requirements/py38_windows_conda.yml

Or use the file `py38_linux_conda.yml` for Linux.


Resolve Dependencies
********************

Dependencies can be installed by letting pip or conda resolve versions from requirements
files. For conda (untested for pip), this can take a long time, so using Mamba_ is recommended.

If you would like a specific PyTorch version to use with PyHDX (ie CUDA/ROCm support), you should install this first.
Installation instructions are on the Pytorch_ website.

Then, install the other base dependencies and optional extras.

To install all dependencies, including development tools:

.. code-block:: rst

    $ conda install --file requirements/req-all.txt

Or choose which extras to install by using the 'req-<extra>.txt' files.

Install PyHDX in develop/editable mode

.. code-block:: rst

    $ conda develop .

Or

.. code-block:: rst

    $ pip install -e


Running from source
-------------------

To run the web application:

.. code-block::

    $ python pyhdx/web/serve.py

This runs the pyhx web application without a Dask cluster to submit jobs to, so
submitting a fitting job will give an error.

To start a dask cluster separately, open another terminal tab and run:

.. code-block:: rst

    $ python pyhdx/local_cluster.py


Configuration
-------------

A configuration file is located in the `.pyhdx` folder in the user home directory. This file
is used by default and can be edited to change PyHDX default settings.

Alternatively, users can create additional `.yaml`configuration files in this directory, after
which the scripts `local_cluster.py` and `serve.py` prompt the user for which file to use.

The section `server` configures the panel server settings. In this section the additional keys
`port` and `websocket_origin` can be added, which are passed to `panel.serve`. See the panel
`Deploy and Export`_ deploy section for more information.


.. _Github repo: https://github.com/Jhsmit/pyhdx

.. _Pytorch: https://pytorch.org/

.. _Mamba: https://mamba.readthedocs.io/en/latest/

.. _Deploy and Export: https://panel.holoviz.org/user_guide/Deploy_and_Export.html
