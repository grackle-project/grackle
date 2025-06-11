Style Formatting
================

The Grackle repository is configured with tools to enforce checks (such as running the linter or applying code formatting) on the various files in the repository.

The primary checks pertain to python code: all python code is linted and all new python files are formatted by automated tools (for the time being, code in older files will not be formatted to avoid merge conflicts).

.. note::

   In the near future (see :gh-pr:`222`), we will start formatting C/C++ code.

Because we realize that these tools may seem overwhelming or complicated, we have configured the repository to try to simplify the experience to the greatest extent possible.
All of these checks are managed by the `pre-commit <https://pre-commit.com/>`__ software.
We will discuss how to automatically invoke these checks as part of a PR (no local installation is required) or on your local machine down below.

Automatic Checks (No installation required)
-------------------------------------------

At a basic level, you don't have to worry about manually invoking any of these tools on your machine.
In fact, you are free to entirely ignore the existence of these tools until it comes time to submit a Pull Request.
When you submit a Pull Request (and whenever you update it), various forms of continuous integration are triggered.

For the present discussion, the `pre-commit.ci <https://pre-commit.ci/>`__ continuous integration tool is of primary relevance.\ [#f1]_
This tool executes all of the formatting tools and if your submission doesn't satisfy all of the requirements, it will fail and report each problem.

If you don't want to install anything locally, there are 2 approaches for addressing problems:

1. You can leave a comment on the Pull Request that simply states

      pre-commit.ci autofix

   and pre-commit.ci will contribute push a commit to your branch that fixes as many issues as possible.
   **This is the recommended way to fix all code-formatting issues.**

2. You can also manually fix the issues locally (and then push your changes).
   **You will need to do this to address most linting errors.**


Running the Checks Locally
--------------------------

To run the checks locally, we strongly encourage you to install the pre-commit software.
This software is written in python and can be installed with ``pip``.
The `installation instructions <https://pre-commit.com/#installation>`__ also mention an alternative approach where you can download and run pre-commit without fully installing it (as a "zipapp").

Once you have installed ``pre-commit``, you can enforce the checks by invoking the following command from the root of your Grackle repository:

.. code-block:: shell-session

   ~/grackle $ pre-commit run --all-files

The above command does 2 things:

1. First, it ensure that local copies of the correct versions of the required enforcement tools are installed.
   These local copies are only accessed by pre-commit and won't affect other parts of your system.
   These copies are also cached (so that the tools don't need to be reinstalled on every invocation).

2. Then the command applies the enforcement tools on the files in your repository (tool-specific exclusions, like files listed by ``tool.ruff.format.exclude`` in **pyproject.toml**, are obviously respected).

.. caution::

   The above command will modify the files in your repository (after all, that's the whole point of the command).
   The pre-commit software does not provide a way to reverse this change.

.. important::

   As already noted, we **strongly** recommend using pre-commit to run the checks locally in order to make sure the version you are using is consistent with the version used in Continuous Integration.
   Manually downloading and invoking a tool like ``ruff``, **probably** won't cause many issues.
   However, using pre-commit will become very important when we start formatting C/C++ code.


Summary of Code Checks
----------------------

Ruff Formatter (Python code formatting)
+++++++++++++++++++++++++++++++++++++++

Python code is formatted by the `Ruff Formatter <https://docs.astral.sh/ruff/formatter/>`__ tool.
This is provided as part of the popular `Ruff <https://github.com/astral-sh/ruff>`__ tool (i.e. it's invoked with ``ruff format``).
It performs a similar role to the `Black code formatter <https://black.readthedocs.io/en/stable/>`__.

Ruff Linter (Python code linting)
+++++++++++++++++++++++++++++++++

Python code is linted by the `Ruff Linter <https://docs.astral.sh/ruff/linter/>`__ tool.
This is also provided as part of the popular `Ruff <https://github.com/astral-sh/ruff>`__ tool (i.e. it's invoked via ``ruff check``).
It performs a similar role to the `Flake8 linter <https://pypi.org/project/flake8/>`__.


Miscellaneous checks
++++++++++++++++++++
Some miscellaneous checks are also implemented by a set of miscellaneous enforcement scripts provided by the authors of pre-commit.

.. rubric:: Footnotes

.. [#f1] It's worth clarifying there are essentially 3 distinct entities named "pre-commit":

   1. the `pre-commit <https://pre-commit.com/>`__ software.
      Contributors **only** need to know about this if they want to apply the enforcement tools locally.
   2. the `pre-commit.ci <https://pre-commit.ci/>`__ continuous integration service.
      This is named because the service simply executes the pre-commit software.
      (All contributors will encounter this)A.
   3. the "pre-commit" `git hook <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>__`.
      This is one of multiple different "hooks" offered by git.
      The pre-commit software is named after this hook because it was originally designed to be used with this hook.
      We do **NOT** currently recommend that contributors use the pre-commit software with the pre-commit hook (unless they fully understand the choice that they are making).
