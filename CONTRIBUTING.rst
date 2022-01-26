.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated!

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/Jhsmit/pyhdx/issues.

If you are reporting a bug, when running PyHDX locally, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug, possibly include a minimal dataset which reproduces the bug.

If you are reporting a bug when using the hosted PyHDX web interface, please include:

* The version of PyHDX as shown in the top header
* Which web component you were using and steps to reproduce the bug.
* Possibly include input/output data which helps reproduce the bug.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/Jhsmit/pyhdx/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.

Pull Requests
~~~~~~~~~~~~~

You can contribute code by submitting a pull request. If you contribute new features describe in your
PR what the new feature does, possible references to literature and how it should be used.

The PR should add new tests for the new feature and all current tests should pass. New functions and
classes should have docstrings and ideally code examples and documentation should be added.

Documentation
~~~~~~~~~~~~~

Docstrings use `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`__ style headers.
Parameters for `__init__` are documented at the `class` level. `intersphinx <http://www.sphinx-doc.org/en/stable/ext/intersphinx.html>`__
is used to link to classes of external packages.
Some guidelines for how to refer to objects or entities:

* Variable, module function and class names: ```numpy```
* Boolean values: ````True````
* Python objects: ``:obj:`int```, ``:obj:`bool``
* Numpy arrays: ``:class:`~numpy.ndarray```
* Pandas dataframe: ``:class:`~pandas.DataFrame```
* Optional array-like:  `` ``
* Default values: ``copy : :obj:`bool`, default: ``True````

Parameters:

.. code::
    Parameters
    ----------
    filename : :obj:`str`
    copy : :obj:`bool`
    dtype : data-type
    iterable : iterable object
    shape : int or tuple of int
    files : list of str

Args and kwargs:

.. code::
    *args : :obj:`tuple`
        Additional arguments should be passed as keyword arguments
    **kwargs : :obj:`dict`, optional
        Extra arguments to `somefunction`.

Properties:

.. code::
    *args : :obj:`tuple`
        Additional arguments should be passed as keyword arguments
    **kwargs : :obj:`dict`, optional
        Extra arguments to `somefunction`.

Future docstring typing will probably be through `PEP 484 <https://www.python.org/dev/peps/pep-0484/>`__
type hints.