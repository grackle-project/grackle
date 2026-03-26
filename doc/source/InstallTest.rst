.. _install-tests:

install-tests
=============

Our install-tests primarily serve as a check on the build-system and on our documentation.
They are driven by the :source:`tests/install-tests/installtest.py` test harness.

At a high-level, a single test-case is quite simple.
The test harness performs four phases:
1. set up a scenario where Grackle has been appropriately set up for the test case (e.g. maybe Grackle has been installed as a shared library)
2. set up a sample project (this might involve copying a **Makefile** or **CmakeLists.txt** file).
3. the test harness tries to execute a sequence of test-case specific shell commands to build a sample program.
4. the harness tries to execute the test program.

Motivation
----------

We (the developers) provide a lot of instruction (e.g. see :doc:`here <Integration>`) and convenience-machinery (e.g. CMakge Package Config Files and pkg-config files) to try to make it as easy as possible to use Grackle in external codes.
In our experience, it is easy to break these instructions or the convenience machinery if it isn't regularly tested.

The reader may understandably wonder "why don't other open-source software library projects have such a big emphasis on these sorts of tests?"
There are a few reasons:

1. There is often an assumption that the only people that will be building a software library are developers who are already well-versed in building/linking/consuming software libraries.
   - When popular libraries (think of hdf5, fftw3, libpng, zlib, xzutils) are used in projects written less experienced developers, the projects are typically made available by package managers [packagers]_ or may distribute header-only libraries (only works if all consumers are written in C++) 
   - This obviously isn't a great assumption for users of astrophysical simulation codes that may depend on Grackle (thus, we make an effort to provide lots of instruction).

2. We attempted to make Grackle easy to use as a static library in pure C programs.
   It turns out that this is at the root of a lot of complexity:
   - complications arise because Grackle itself has external dependencies
     - it explicitly depends on hdf5.
     - it also has implicit dependencies a regular C compiler (or linker) won't know about out of the box:
       - Most obviously, this is the C++ standard library (historically, it was the Fortran standard library).
       - When compiled with OpenMP, Grackle has an additional implicit dependency on the library.
   - The fact that we also took steps to make it easy to install both shared and static copies of Grackle to a single directly adds additional complexity.
   - Moreover, we added slightly more complexity by making it possible for external CMake projects more easily select between the use of shared/static libraries

   .. note::

      To be clear, the only real alternative was to act like static libraries don't exist outside of embedded Grackle builds (within external CMake projects).

3. Grackle doesn't have a stable ABI and we attempt to remain somewhat backwards compatible with choices made by the classice build-system.
   In practice, this translates to creating a ``libgrackle.so`` (or ``libgrackle.dylib`` file on MacOS) symlink during installation that links to a somewhat atypically named shared library during the install step.
   While the logic for doing this is simple, it should be tested.


Why not run these tests as part of the core test suite?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At a surface level, the steps of a single test-case makes it sound like these tests could be part of the core-library test suite.
In fact, the core-lib tests already involve a set of "compile-tests" that closely resemble the structure of some tests.
Moreover, other open source projects might try to test this sort of thing by constructing a ctest test-case that explicitly calls ``ctest`` executable with the ``--build-and-test`` option (an example can be found `here <https://github.com/fmtlib/fmt/blob/dc05bee30755bb993add7ed90845d87c2315f9c6/test/CMakeLists.txt#L166>`__).
There are 2 limitation to this approach:
1. it **only** tests the direct product of a build (if you build Grackle as a shared library, you only test builds against a shared library)
2. tests commonly try to link against the target that's in the build-directory (when in practice people may want to link against Grackle after it has been fully installed).\ [link-test]_

These limitations aren't usually as much of a concern in external project because (as noted above) there is usually far less emphasis on making the project easy to install/link against, for less experienced developers.


Running the tests
-----------------

Running the tests is quite simple:

1. Make sure ``docker`` is installed.
   On MacOS, you can use docker-desktop.
   An explanation for why we use ``docker`` is provided :ref:`below <installtest-why-docker>`.

2. Simply execute the :source:`tests/install-tests/installtest.py` script.

.. note::

   If you are using docker-desktop, you may need to actually launch the Docker Desktop GUI, before you invoke the ``installtest.py`` test harness (the Docker Desktop GUI needs to keep running for the duration of the tests).
   If you forget to do this, you may get cryptic error messages.

   It's a little unclear at this time whether this is always a requirement or only a requirement when docker-desktop is installed in a particular way (maybe it's only a thing on MacOS and not on Linux -- after all Docker Desktop on MacOS needs to spin up a Linux VM behind the scenes).


All of these steps are completed within a ``docker container``.


Design Rationale
----------------


Big Picture Ideas
+++++++++++++++++

While designing the test harness, we had the following (somewhat related) guiding principles:

- the tests should not require gracklepy to be installed
  * Ideologically, these tests pertain to the core library and have nothing to do with ``gracklepy``
  * Practically, there's value to being able to run the tests without installing ``gracklepy``.
    Specifically, these tests explicitly check for the exact kinds of issues that would break builds of ``gracklepy``.

- these tests should be portable and easy to run
  * practically, this means that the runner should use the system provided python installation
  * currently, we the tests just require the installation of docker (while that's a bit of a hurdle, it's :ref:`somewhat unavoidable <installtest-why-docker>`)
  * while I'm not opposed to using external python packages in the test harness, we would probably want to retain a separate pyproject.toml file (or at least a REQUIREMENTS.txt file) for the harness.

- retain the capacity to easily integrate these tests with CTest
  * if the tests were "integrated," then CTest could be used execute individual or all of the tests.
    This could make the process of launching tests related to the core library more convenient.
  * for added clarity, all of our googletest unit tests are integrated with CTest (for that precise reason, the test harness's CLI takes cues from googletest's interface).
  * at present, we have not performed this integration, but we still want to retain the capacity to do it in the future.\ [ctest-integration]_

Why don't we use pytest?
++++++++++++++++++++++++

The short answer is: it wouldn't really help us.

To elaborate:
- pytest isn't really designed for these kinds of tests; it's primarily designed to test python code.
- Since pytest is extendible, we could definitely extend it to make it work for this purpose.
  In practice, we would have to write much of the code we already have.
  Using pytest might realistically save us 100-200 lines of codes.
- But this comes with a few drawbacks:
  * this choice sacrifices the ability to easily integrate the install-tests with CTest (this may not a big deal, but I'm a little hesitant to make that choice).
  * Refactoring the creation of docker images to use the fixture system would take some effort (it's definitely doable).
  * For non-experts, pytest's control flow is fairly tricky to follow.
  * Least importantly, as the installtest.py's first external python dependency, using it would require us to give some thought about specifying python dependencies (not a big deal)

Why we use TOML files for configuring tests
+++++++++++++++++++++++++++++++++++++++++++

We decided:
- It *probably* makes more sense to enumerate the testing cases with TOML files than describing each and every test in python code.
- It's important that we be able to embed explanatory comments into the config files
  (this rules out JSON files).

This leaves YAML, INI, and TOML files.
I feel strongly that INI is the wrong choice since its semantics are **NOT** well-defined.

TOML is preferable to YAML since TOML is "just" sophisticated enough to meet our needs.
While YAML is significantly more powerful, it is also significantly more complex (it's extremely easy to make subtle mistakes).

TOML is also preferable from a dependency perspective:
- Currently, we include :source:`tests/install-tests/vendored_tomli.py` which is a vendored version of the `tomli <https://pypi.org/project/tomli/>`__ package (the file is ~800 lines).
- Although vendoring is obviously undesirable, it's very much a short term solution.
  Once we are willing to assume everybody has a minimum python version of 3.11, we can just use python's builtin ``tomllib`` package.
- In comparison, vendoring ``PyYAML`` is simply out of the question


.. _installtest-why-docker:

Why use Docker?
+++++++++++++++


.. note::

   ADD ME



How the tests work
------------------

In practice, step 1 consists of creating a docker image, while steps 2-4 are executed in a docker container based on the image.
The use of docker is intended to take adva

.. note::

   ADD ME

Adding a new test case
----------------------

.. note::

   ADD ME


.. rubric:: Footnotes

.. [packagers] In particularly popular projects that are very wildly used by less experienced developers , a group of non-developers known as packagers (e.g. Debian or Homebrew packagers) take charge of distributing pre-compiled copies of the software.

.. [link-test] It *is* possible to temporarily overwrite the install-directory to try to temporarily install Grackle to that directory (by invoking ``cmake --install --prefix ``<DEST>``) and then test linking against the contents of that directory.
   However there are 2 wrinkles to that approach.
   First, not all install-tests will work very well in this case (think about install-tests that try to link a test-program against Grackle without using CMake or pkg-config).
   Second (and more importantly), different CMake logic is executed when the install destination is specified at initial configuration vs using the temporary overwrite logic (we have in fact had bugs pertaining to this logic).

.. [ctest-integration] To integrate with CTest, we would need to make sure that we properly label each installtest test case and possibly disable them when a simple ``ctest`` command is invoked since these tests take a lot longer than the other kinds of tests (this is doable, but would take a little effort).
   We would also need to properly record each test case's dependency on the task of building a docker image (again, very doable).
