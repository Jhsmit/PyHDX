"""Module of Invoke test tasks for testing the CODE QUALITY to be invoked from the command line. Try

invoke --list=test

from the command line for a list of all available commands.
"""

from invoke import task


@task()
def bandit(command):
    """Runs Bandit the security linter from PyCQA."""
    print(
        """
Running Bandit the Python Security Linter
to identify common security issues in Python code
=================================================
"""
    )
    command.run("bandit -r ./", echo=True)


@task()
def black(command):
    """Runs black (autoformatter) on all .py files recursively"""
    print(
        """
Running Black the Python code formatter
=======================================
"""
    )
    command.run("black .", echo=True)


@task()
def flake8(command):
    """Runs flake8 linter."""
    print(
        """
Running flake8.
Flake8 is a python tool that glues together pycodestyle, pyflakes, 
mccabe, and third-party plugins to check the style and quality of 
some python code. 
=======================================================================
"""
    )
    command.run("flake8 .", echo=True)


@task()
def isort(command):
    """Runs isort (import sorter) on all .py files recursively"""
    print(
        """
Running isort the Python code import sorter
===========================================
"""
    )
    command.run("isort .", echo=True)


@task
def pytest(
    command,
    test_files="tests",
):
    """Runs pytest to identify failing tests

    Arguments:
        command {[type]} -- Invoke command object

    Keyword Arguments:
        test_files {str} -- A space separated list of folders and files to test. (default: {'tests})


    # Print running pytest
    """
    print(
        """
Running pytest the test framework
=================================
"""
    )
    # Build the command_string
    command_string = f"pytest {test_files}"

    # Run the command_string
    command.run(command_string, echo=True)



@task()
def pylint(command, files="setup.py tasks pyhdx tests"):
    """Runs pylint (linter) on all .py files recursively to identify coding errors

    Arguments:
        command {[type]} -- [description]
        files {string} -- A space separated list of files and folders to lint
    """
    # https://stackoverflow.com/questions/22241435/pylint-discard-cached-file-state
    # from astroid import MANAGER
    # MANAGER.astroid_cache.clear()
    print(
        """
Running pylint.
Pylint looks for programming errors, helps enforcing a coding standard,
sniffs for code smells and offers simple refactoring suggestions.
=======================================================================
"""
    )
    command_string = f"pylint {files}"
    command.run(command_string, echo=True)


@task
def mypy(command, files="setup.py tasks pyhdx tests"):
    """Runs mypy (static type checker) on all .py files recursively

    Arguments:
        command {[type]} -- [description]
        files {string} -- A space separated list of files and folders to lint
    """
    print(
        """
Running mypy for identifying python type errors
===============================================
"""
    )
    command_string = f"mypy {files}"
    command.run(command_string, echo=True)


@task(
    pre=[isort, black, bandit, flake8, pylint, mypy, pytest],
    aliases=["pre_commit", "test"],
    name="all",
)
def _all(command):  # pylint: disable=unused-argument
    """Runs isort, autoflake, black, pylint, mypy and pytest

    Arguments:
        command {[type]} -- [description]
    """
    # If we get to this point all tests listed in 'pre' have passed
    # unless we have run the task with the --warn flag
    if not command.config.run.warn:
        print(
            """
All Tests Passed Successfully
=============================
"""
        )
