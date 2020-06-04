.. highlight:: shell

============
Installation
============


Stable release
--------------

(Currently no stable release available. This section will updated soon)

To install PyHDX, run this command in your terminal:

.. code-block:: console

    $ pip install pyhdx

This is the preferred method to install PyHDX, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for PyHDX can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/Jhsmit/pyhdx

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/Jhsmit/pyhdx/tarball/master

pyHDX can then be installed with ``conda`` (requires ``conda build``):

.. code-block:: console

    $ conda develop pyhdx

or ``pip``:

.. code-block:: console

    $ pip install pyhdx


..
    Once you have a copy of the source, you can install it with:

    .. code-block:: console

        $ python setup.py install


To launch the web application:

.. code-block::
    $ panel serve panel/main.py

.. _Github repo: https://github.com/Jhsmit/pyhdx
.. _tarball: https://github.com/Jhsmit/pyhdx/tarball/master
