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

test

#. Step 1.
#. Step 2.
#. Step 3.

testtest

* This is a bulleted list.
* It has two items, the second
  item uses two lines.

testset

-  This is a bulleted list.
-  It has two items, the second
   item uses two lines.

PyHDX web application can be launched from the command line using ``pyhdx`` command with below options,

* To run PyHDX server using default settings on your local server,

        ``pyhdx serve``

* To run PyHDX server using the IP address and port number of your dask cluster,

        ``pyhdx --cluster <ip> <port>``


Dependencies
------------

The requirements for PyHDX are listed in requirements.txt

.. _Github repo: https://github.com/Jhsmit/pyhdx
