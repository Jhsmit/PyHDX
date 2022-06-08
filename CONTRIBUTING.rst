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

Docstrings use `google <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`__ style headers.
Parameters for `__init__` are documented at the `class` level. `intersphinx <http://www.sphinx-doc.org/en/stable/ext/intersphinx.html>`__
is used to link to classes of external packages.
Some guidelines for how to refer to objects or entities:

* Variable, module function and class names: ```numpy```
* Boolean values: ````True````
* Python objects: ``:obj:`int```, ``:obj:`bool``
* Numpy arrays: ``:class:`~numpy.ndarray```
* Pandas dataframe: ``:class:`~pandas.DataFrame```
* Classes within the same module: ``:class:.HDXMeasurement``
* Default values: ``copy : :obj:`bool`, default: ``True````
* Refer to arguments / parameters with single ticks: ```param1```


.. code::
    class Coverage(object):
    """Single docstring line description.

    Args:
        data: Dataframe with input peptides.
        c_term: Residue index number of the C-terminal residue (where first residue in index number 1).
        * args: Additional arguments
        **metadata: kwargs named metadata.


    """

pahtlike should be: Pathlike[str] or Union[str, Pathlike[str]]

Properties:

.. code::
    *args : :obj:`tuple`
        Additional arguments should be passed as keyword arguments
    **kwargs : :obj:`dict`, optional
        Extra arguments to `somefunction`.

Additional resources:

`Sphinx docstrings <https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html>`__
`Mypy cheatsheet <https://mypy.readthedocs.io/en/latest/cheat_sheet_py3.html>__
`Google styleguild <https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings>__