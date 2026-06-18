Formatting, Linting, Other Checks
=================================

The Grackle repository is configured with tools to enforce checks (such as running the linter or applying code formatting) on the various files in the repository.

The primary checks pertain to python and C/C++ code.
All code is linted and all new code files are formatted by code-formatters.
For the time being, code in older files will not be formatted to avoid merge conflicts.

Because we realize that these tools may seem overwhelming or complicated, we have tried to simplify the process of using them to the greatest extent possible:

1. All checks are automatically run (through continuous integration) when you make or modify a pull request.
   Furthermore, we support machinery that contributors can trigger from the GitHub Pull Request web interface to add commits their your code in order to satisfy failed code-formatting checks (manual intervention is generally required to address failed linting checks).
   This **makes it possible to use these tools without locally installing them.**
   We provide more information :ref:` down below <pr-auto-checks>`

2. Most of the checks are managed by the `pre-commit <https://pre-commit.com/>`__ software.
   All such check are very easy to install and run locally (details are provided :ref:`below <local-pre-commit-checks>`).
   We also talk a little about running the other more-manual checks below.
   
At the bottom of this page, we provide :ref:`an overview of all checks <other-static-code-analysis>` that the Grackle project uses.

.. _pr-auto-checks:

Automatic Checks (No installation required)
-------------------------------------------

At a basic level, you don't have to worry about manually invoking any of these tools on your machine.
In fact, you are free to entirely ignore the existence of these tools until it comes time to submit a Pull Request.
When you submit a Pull Request (and whenever you update it), various forms of continuous integration are triggered.

The `pre-commit.ci <https://pre-commit.ci/>`__ continuous integration service is responsible for triggering most of the checks (i.e. the same checks that we would :ref:`run locally <local-pre-commit-checks>` with the pre-commit software).\ [#f1]_
At the time of writing, the remaining check(s) are run by workflow(s) the CircleCI continuous integration services.

These tools execute all of the checks and if your Pull Requests submission doesn't satisfy all checks the relevant continuous integration workflow(s)/task(s) will report failure(s).

If you don't want to install anything locally, there are 2 approaches for addressing problems:

1. You can leave a comment on the Pull Request that simply states

      pre-commit.ci autofix

   and pre-commit.ci will push a commit to your branch that fixes as many issues as possible.
   **This is the recommended way to fix all code-formatting issues.**

2. You can also manually fix the issues locally (and then push your changes).
   **You will need to do this to address most linting errors.**

.. _local-pre-commit-checks:

Locally run checks registered with Pre-Commit
---------------------------------------------

If you are interested in running one or more check locally, we **STRONGLY** encourage you to install the pre-commit software and use pre-commit to invoke checks rather than manually invoking the underlying tool, directly.
(We provide separate instructions :ref:`in the next section <other-static-code-analysis>` for tools that aren't registered with pre-commit)

.. important::

   It is particularly important that you invoke ``clang-format`` through the pre-commit software.
   Please read our :ref:`section on clang-format <clang-format>` for more details.

The full list of checks (and some configuration options) can be found in the :source:`.pre-commit-config.yaml` file.

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


.. _other-static-code-analysis:

Locally run other Static Code Analysis Checks
---------------------------------------------

For each tool not registered with the pre-commit framework, we rely upon built-in CMake support since the tool must be passed a subset of the compiler flags, have access to an auto-generated header file, or both.
Currently, to locally run such tools, you need to (unfortunately) take some manual steps.

First, we'll cover the generic procedure to use such tools (fortunately, it's very generalizable), and then we'll list the actual tools in this category that we actually use.

Generic Procedure to Run Tools
++++++++++++++++++++++++++++++

CMake adopts a co-compilation approach to support these kinds of static analysis tools.
Essentially, we configure an ordinary\ [#f2]_ build of the core Grackle library (examples and unit-tests should generally be enabled) and specify an extra configuration option for each tool of interest (details given :ref:`below <summary-of-checks>` on a tool-by-tool basis).
Then, we execute a normal build.
When the build-system executes the compiler to compile a C++ source file into an object file, the build-system launches each of the registered tools, and any failing checks are reported as if they are compiler errors.
If there aren't errors, CMake will successfully build all relevant targets.

The preceding procedure is sketched out below for an arbitrary build directory, ``<build-dir>``:

.. code-block:: shell-session

   ~/grackle $ cmake <normal-args...> <tool-args...> -B <build-dir>
   ~/grackle $ cmake --build <build-dir>

CMake's support for these tools implicitly assumes that these tools are relatively fast.
When that assumption is true, the premise is that you would configure all relevant tools in a build-directory that you are using during development, and everytime you recompile after making a change, all of the tools would also run.
In practice, it's not necessarily true (e.g. as clang-tidy checks more and more things, it takes longer and longer).
Consequently, we recommend that you make a separate build-directory for running the checks.

.. note::

   Because these tools are run alongside the compiler and CMake supports incremental rebuilds, there are some noteworthy consequences.
   Within a given build directory, `<build-dir>`:

   * once a source file has been succesfully compiled, it will not be compiled again (the static analysis tool won't be run again) until after the source file (or an included header) has been modified.
   * at this time, modifying a static analysis tool's configuration file won't force the tool to be run again
   * to force the tool to run again, you can always build the clean target first. 
     This can be done by invoking `make clean` or `ninja clean` within the build directory (depending on whether you use a Makefile or Ninja generator).
     Or, you could use `cmake --build <path/to/build-dir> --target 'clean'`


Tools in this category
++++++++++++++++++++++

We provide a list of each tool in this category (with links to further details) down below:

1. clang-tidy (:ref:`more details <clang-tidy>`)

Yes, there's currently just 1 tool, but we could add support for others in the future (like `include-what-you-use <https://include-what-you-use.org/>`__)


.. _summary-of-checks:

Summary of All Code Checks
--------------------------

.. _clang-format:

``clang-format`` (C/C++ Formatting)
+++++++++++++++++++++++++++++++++++
C/C++ code is formatted by `clang-format <https://releases.llvm.org/18.1.8/tools/clang/docs/ClangFormat.html>`__.

.. important::

   The ``clang-format`` version number is tied to the version number of the entire LLVM project, which has a fairly rapid release cadence and clang-format is **NOT** forward or backwards compatible.
   If you want to manually install and invoke clang-format on your machine (outside of the pre-commit framework), you **NEED** to make sure that you use the exact same version of clang-format that is used by the pre-commit framework.
   If you use a different version, differences **will** arise,\ [#f3]_ (we may ask you to rewrite your Pull Request's commit-history to remove lots of spurious commits that change this file due to this mismatch).

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


.. _clang-tidy:

``clang-tidy`` (C++ Linting)
++++++++++++++++++++++++++++

C++ code is linted by `clang-tidy <https://clang.llvm.org/extra/clang-tidy/>`__.

Details about the enforced checks are specified in the :source:`.clang-tidy` file.
We go into more detail about executing ``clang-tidy``, adjusting the list of checks, and getting an appropriate version of ``clang-tidy`` down below.

Executing the tool
^^^^^^^^^^^^^^^^^^

As noted :ref:`above <other-static-code-analysis>`, ``clang-tidy`` is *NOT* registered as part of pre-commit framework, since we rely upon CMake's built-in support for running the tool.
The linked section also provides general guidelines for using a tool like ``clang-tidy``.

We support the ``GRACKLE_DEV_CLANG_TIDY`` cmake configuration variable as a convenient way to configure clang-tidy.
This should specify a semicolon-delimited list.
The first list element specifies the name of the ``clang-tidy`` executable (CMake searches for the executable in all the standard places) or the absolute path to the executable.
Any subsequent list elements specify command-line arguments are forwarded to the executable.
For example, we might pass ``-D 'GRACKLE_DEV_CLANG_TIDY=clang-tidy;-warnings-as-errors=*'`` to indicate that clang-tidy should report warnings as errors.
This variable works around a bug that can come up when using ``clang-tidy`` (and actually backports some logic for CMake versions older than 3.25 to make this work).

More experienced CMake users are free to directly specify the standard `CMAKE_<LANG>_CLANG_TIDY <https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_CLANG_TIDY.html>`__ variable.

What Happens When I Encounter a Spurious Failure?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``clang-tidy`` documentation `has a section <https://clang.llvm.org/extra/clang-tidy/#suppressing-undesired-diagnostics>`__ on suppresing undesired diagnostics.
You can add comments to disable clang-format on sections of code (if you do that, you should probably add an explanatory comment).

``clang-tidy`` has a lot of checks and some of them can be very opinionated.
If a particular check feels particularly burdensome or "keeps getting in the way", you should feel free to disable it.

This is especially true at the time of writing because
- the initial configuration is very liberal about enabling checks (it's very plausible that we unintentionally enabled something overly restrictive).
- we are still in the process of converting the codebase to C++ (the priority is get working C++ code first, and then refactor after -- we shouldn't let checks slow down this process)

.. important::

   Before disabling a check, please confirm that the :source:`.clang-tidy` file has doesn't have a comment (usually after the list of checks) providing a strong reason why that check should be enabled.

Getting an appropriate version of clang-tidy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While ``clang-tidy`` generally has a reputation for being more stable between versions than ``clang-format`` it's probably advisable to try to match the major version of ``clang-tidy`` that is currently used in continuous integration.

There are a few ways to get a working copy of clang-tidy with a given ``<MAJOR_VERSION>`` (like ``19``, ``20``, etc.) without building explicitly building from source (which is actually a pretty seamless experience):

.. tabs::

   .. tab:: use ``apt`` (Debian/Ubuntu)

      On Debian or Ubuntu, you can use the `LLVM provided packages <https://apt.llvm.org/>`__.
      After adding the archive signature to ``apt`` (instructions provided on the linked website), you would invoke

      .. code-block:: shell-session

         $ sudo apt-get install clang-tidy-<MAJOR_VERSION>

      At the time of writing, this places an executable named ``clang-tidy-${MAJOR_VERSION}`` into your path

   .. tab:: use homebrew (macOS)

      On macOS (and linux), you could use the homebrew formula to install clang-tidy as part of the llvm package.
      At the time of writing, after you invoke

      .. code-block:: shell-session

         $ brew install llvm@<MAJOR_VERSION>``

      you can find the relevant executable at ``<hbrew-prefix>/opt/llvm@<MAJOR_VERSION>/bin/clang-tidy`` (the ``<hbrew-prefix>`` is ``/opt/homebrew`` for Apple Silicon and ``/usr/local`` for macOS Intel)

  .. tab:: via PyPI (most general)

     On any Unix-like OS you can get through ``clang-tidy`` `from PyPI <https://pypi.org/project/clang-tidy/>`__.
     On main-stream operating systems and architectures, you can download a precompiled copy (other platforms, your python package-manager will try to compile the executable).

     To install ``clang-tidy`` into a virtual environment at ``<path/to/venv>``, you might invoke (there are multiple ways to do this) the following commands:

     .. code-block:: shell-session

        $ python3 -m venv <path/to/venv>
        $ <path/to/venv>/bin/python3 -m pip install 'clang-tidy==<MAJOR_VERSION>.*'

     At the time of writing, these commands install the relevant executable at ``<path/to/venv>/bin/clang-tidy``.
     In this scenario, uninstalling the tool is as simple as using ``pip uninstall`` or deleting the venv.


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

.. [#f2] Technically, this only works in CMake builds that use a Makefile or Ninja generator.
      In just about any context relevant to Grackle, those are the only relevant generators.

.. [#f3] While this hasn't been tested, it's *probably* okay to use a different version of clang-format as long as the Major and Minor version numbers match (e.g. if the pre-commit plugin was configured to use version 20.1.2 of clang-format, you could probably use version 20.1.6).
