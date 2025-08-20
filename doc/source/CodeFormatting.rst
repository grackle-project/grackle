Style Formatting
================

The Grackle repository is configured with tools to enforce checks (such as running the linter or applying code formatting) on the various files in the repository.

The primary checks pertain to python and C/C++ code.
All python code is linted and all new python files are formatted by automated tools.
Likewise all new C/C++ code in new files are formatted by automated tools.
For the time being, code in older files will not be formatted to avoid merge conflicts.

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

To run the checks locally, we **STRONGLY** encourage you to install the pre-commit software and use pre-commit to invoke the checks, rather than manually invoking linters or formatters directly.

.. important::

   It is particularly important that you invoke ``clang-format`` through the pre-commit software.
   Please read our :ref:`section on clang-format <clang-format>` for more details.

The pre-commit software is written in python and can be installed with ``pip`` .
The `installation instructions <https://pre-commit.com/#installation>`__ also mention an alternative approach where you can run download and run pre-commit without fully installing it (as a "zipapp").

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


Summary of Code Checks
----------------------

.. _clang-format:

``clang-format`` (C/C++ Formatting)
+++++++++++++++++++++++++++++++++++
C/C++ code is formatted by `clang-format <https://releases.llvm.org/18.1.8/tools/clang/docs/ClangFormat.html>`__.

.. important::

   The ``clang-format`` version number is tied to the version number of the entire LLVM project, which has a fairly rapid release cadence and clang-format is **NOT** forward or backwards compatible.
   If you want to manually install and invoke clang-format on your machine (outside of the pre-commit framework), you **NEED** to make sure that you use the exact same version of clang-format that is used by the pre-commit framework.
   If you use a different version, differences **will** arise,\ [#f2]_ (we may ask you to rewrite your Pull Request's commit-history to remove lots of spurious commits that change this file due to this mismatch).

   This version number is stored in the ``.pre-commit-config.yaml`` file by the ``rev`` parameter in the section for the "clang-format plugin".

   For this reason, we **STRONGLY** encourage you to invoke the ``pre-commit`` tool for local formatting (the pre-commit tool will install a sandboxed version of ``clang-format`` that has the correct version).

About ``clang-format``:

* At the time of writing, ``clang-format`` enforces formatting rules that are largely derived from the google style (with a handful of tweaks that derive from the llvm style guide).
  Details about the enforced style are configured in the ``.clang-format`` file at the root of the Grackle repository.

* files that are formatted this way will generally have far fewer merge conflicts.

* Adherence to the formatting guidelines enforced by ``clang-format`` is a requirement for acceptance of pull requests, and deviation from that will require you to explicitly justify why your code is not adhering to that format.
  As already noted, we have temporarily disabled clang-format from applying to most older files (listed in ``.clang-format-ignore``), that existed before we adopted clang-format in order to minimize merge conflicts (we plan to enable ``clang-format`` for these files in the future).

* **NOTE:** Trying to manually learn all of the style-rules is an exercise in frustration.
  Instead, we recommend that you rely upon clang-format.

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

.. [#f2] While this hasn't been tested, it's *probably* okay to use a different version of clang-format as long as the Major and Minor version numbers match (e.g. if the pre-commit plugin was configured to use version 20.1.2 of clang-format, you could probably use version 20.1.6).
