.. highlight:: shell

============
Installation
============


Stable release
--------------

(Currently PyHDX is not yet distributed on PyPi/conda. This section will updated soon)




From sources
------------

The sources for PyHDX can be downloaded from the `Github repo`_.

To install the dependencies:

.. code-block:: console

    $ conda install --file requirements-conda.txt -c conda-forge -c pyviz -c pytorch

Then clone the repository:

.. code-block:: console

    $ git clone git://github.com/Jhsmit/pyhdx

And install pyHDX with ``conda`` (requires ``conda build``):

.. code-block:: console

    $ conda develop pyhdx

To launch the web application:

.. code-block:: console

    $ panel serve panel/serve.py


Dependencies
------------

The requirements for PyHDX are listed in requirements-conda.txt

.. _Github repo: https://github.com/Jhsmit/pyhdx
