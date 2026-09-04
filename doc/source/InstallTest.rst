.. _installtest:

installtest.py
==============

Our install-tests primarily serve as a check on the build-system and on our documentation.
They are driven by the :source:`tests/install-tests/installtest.py` test program.

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


Why don't other open-source projects place emphasis on these tests?
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The reader may understandably wonder "why don't other open-source software projects have such a big emphasis on these sorts of tests?"

To start, it doesn't make sense to compare Grackle with projects that primarily ship executables (e.g. ``git``).
It simply doesn't make sense to consider linking against an executable.

Let's turn our attention to other projects that ship libraries.
Compared to such projects, we still place far more emphasis on these kinds of tests.
This is for several reasons:

1. Developers of other open-source software libraries typically don't expect non-experts to compile their projects.

   .. COMMENT
      We touch on this point (i.e. why we care more about describing the installation process and making it as easy as possible) more than other projects in :ref:`cmake_buildsystem_design_rationale`.
      If we feel compelled to touch on this point again, I really think we should relocate this explanation, and just link to it.

   * Developers commonly assume that the primary people that compile their software library are already well-versed in the practice of building/linking/consuming software libraries.
     These people are typically other contributors to the project or developers of a downstream application that consumes the library.
     Particularly popular libraries (think of hdf5, fftw3, libpng, zlib, xzutils) may also be compiled by package maintainers in order to make the software available through package managers (e.g. apt/dpkg, yum/rpm, brew).

   * Under this system, the package maintainers make it possible for less experienced users to get libraries without needing to compile the library without compiling it themselves.
     At this time, no maintainer has expressed any interest in packaging Grackle for us.\ [#packagers]_

   * In contrast, we assume that users of downstream codes probably need to manually install Grackle themselves, and they lack a lot of experience pertaining to the process of compiling/linking software.
     In fact, Grackle could be the very piece of software an undergrad physics/astronomy major or first-year grad student could compile.
     Thus, we make an effort to provide lots of instruction.

2. When a developer **does** want to make a project easy to distribute, they commonly distribute their project either:

   1. as a header-only library OR

   2. as a single source and header library

   Both cases allow a library's source code to embedded within the repository of a downstream application.\ [#headeronly-lib]_
   The issues checked by each installtest aren't very relevant for libraries consumed in one of these ways.
   Unfortunately this strategy has several shortcomings:

   * There is a strong temptation for downstream developers to fix the bugs and make improvements to the copy of the code vendored in downstream application's code repository (without contributing them back upstream).

   * This strategy only works for downstream application written in a compatible version of the of the language that was used to write the library.
     In other words, if we distributed Grackle in this way, it would only be beneficial to simulation codes written in the same version of C++ or newer (at the time of writing, Grackle is written in C++ 17).

   * This strategy doesn't work very well when the library has external dependencies because the dependencies need to be properly communicated with the build system of the downstream application.
     While this challenge is surmountable with hdf5 or openmp, it won't really work be viable as we start introducing GPU support.

   * Additionally, Grackle's need for data files introduces additional complexity.

3. We attempted to make Grackle easy to use as a static library in pure C programs.
   It turns out that this choice produced a lot of complexity:

   * complications arise because Grackle itself has external dependencies

     - it explicitly depends on hdf5.

     - it also has implicit dependencies a regular C compiler (or linker) won't know about out of the box:

       * Most obviously, this is the C++ standard library (historically, it was the Fortran standard library).

       * When compiled with OpenMP, Grackle has an additional implicit dependency on the library.

   * The fact that we also took steps to make it easy to install both shared and static copies of Grackle to a single directly adds additional complexity.

   * Moreover, we added slightly more complexity by making it possible for external CMake projects to more easily select between the use of shared/static libraries.

   .. note::

      To be clear, the only real alternative to this choice of supporting Grackle as a static libraries is tell people that they are on their own outside of embedded Grackle builds (within external CMake projects).

4. Grackle doesn't have a stable ABI and we attempt to remain somewhat backwards compatible with choices made by the classic build-system.
   In practice, this translates to creating a ``libgrackle.so`` (or ``libgrackle.dylib`` file on MacOS) symlink during installation that links to a somewhat atypically named shared library during the install step (the source shared library is named to communicate libgrackle ABI incompatability).
   While the logic for doing this is simple, it needs testing (it has held bugs before).


Why not run these tests as part of the core test suite?
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

At a surface level, the steps of a single test-case makes it sound like these tests could be part of the core-library test suite.
In fact, the core-lib tests already involve a set of "compile-tests" that closely resemble the structure of some installtest test-cases.
Moreover, other open source projects might try to test this sort of thing by constructing a ctest test-case that explicitly calls the ``ctest`` executable with the ``--build-and-test`` option (an example can be found `here <https://github.com/fmtlib/fmt/blob/dc05bee30755bb993add7ed90845d87c2315f9c6/test/CMakeLists.txt#L166>`__).

There are several limitation to trying to move the tests core-lib test suites.

1. The core-lib test suite **only** tests the direct product of a build (if you build Grackle as a shared library, you only test builds against a shared library).
   If we moved these tests into the core-lib test suite, then a single invokation of the core-lib suite would would only be able to run a subset of the cases covered by the existing installetest test cases.
   To run a different subset, you would need to completely recompile Grackle before launching relaunching the core-lib suite.

2. CMake/CTest is designed to drive tests involving build-artifacts that are located in the build-directory.
   Because most people link against Grackle after it has been fully installed, some of the installtest test cases involve consuming Grackle after it has been fully installed.
   There wouldn't be a good way to make the core-lib test suite directly run these tests.\ [#link-test]_

These limitations aren't usually as much of a concern in external project because (as noted above) there is usually far less emphasis on making the project easy to install/link against, for less experienced developers.

.. _installtest-execution:

Running the tests
-----------------

Running the tests is quite simple:

1. Make sure ``docker`` is installed.
   On MacOS, you can use docker-desktop.
   An explanation for why we use ``docker`` is provided :ref:`below <installtest-why-docker>`.

2. Simply execute the ``exec`` subcommand of :source:`tests/install-tests/installtest.py` script:

   .. code-block:: shell-session

      python3 <path/to/installtest.py> exec

.. note::

   If you are using docker-desktop, you may need to actually launch the Docker Desktop GUI, before you invoke the ``installtest.py`` test harness (the Docker Desktop GUI needs to keep running for the duration of the tests).
   If you forget to do this, you may get cryptic error messages.

   It's a little unclear at this time whether this is always a requirement or only a requirement when docker-desktop is installed in a particular way (maybe it's only a thing on MacOS and not on Linux -- after all Docker Desktop on MacOS needs to spin up a Linux VM behind the scenes).


All of these steps are completed within a ``docker container``.

.. _how-installtest-works:

How installtest.py works
------------------------

At a high-level, you can think of :source:`tests/install-tests/installtest.py`, as a highly simplified version of a program for running CI logic that can be executed locally.
Like the GitHub Actions runner (or the CircleCI runner), it sets up prestine testing environments, infers commands to run from config files (we explain why TOML over YAML :ref:`here <installtest-why-toml>`), and reports failed commands.

In more detail, when you actually execute :source:`tests/install-tests/installtest.py` (in ``exec`` mode), there are 2 concrete phases: (i) task collection and (ii) task execution.

Task Collection
+++++++++++++++

Immediately after you start :source:`tests/install-tests/installtest.py`, the program starts scanning the ``scan-dir`` (by default, this is the :repository-dir:`tests/install-tests/cases` directory).
The program looks for all subdirectories in the ``scan-dir``, and it determines the suite of test cases associated with each subdirectory by reading the contained ``conf.toml`` file.

We will discuss exactly how tests are defined :ref:`later <installtest-defining-tests>`.
For now, it's important to understand that each test-case is run in a docker container.
During collection, we record the names of each required docker image.

After the scanning the entirety of the ``scan-dir``, a task list is constructed that consists of building each required docker image and running each test.


Test Execution
++++++++++++++

After all tasks are collected, :source:`tests/install-tests/installtest.py` starts executing each task one-by-one.

When it comes time to execute a test, the program:

1. Creates a new docker container from a previously built docker image.
   The container is configured such that there is a user named ``gr-user``.
   All of the available images hold the Grackle repository (its a subdirectory of the home-directory, named **grackle**).

2. Copies of relevant files to create a source directory for a sample project (the sample project is a subdirectory of the home-directory named **sample-project**).

3. The program executes each specified build command for the test.
   Effectively, it instructs the shell within the container to execute each command, one-by-one.
   If any command should fail, the test is ended.
   This early exit is normally a failure, unless a test explicity expects a command to fail.

4. The program checks that the sample program has successfully been built and can be executed.

Once the test ends, the program cleans up: it deletes the container created for running the test.

.. _installtest-docker-images:

Docker Images
-------------

As we've noted in multiple other places, each individual install-test is run inside of a Docker Container.
We provide further justification for this choice :ref:`here <installtest-why-docker>`.

We identify each image by the name of the build step within :source:`tests/install-tests/installtest.Dockerfile` that the image is built from.
We list each image down below:


.. embed-cli-output:: ../../tests/install-tests/installtest.py
   :args: entrydoc image-target

.. _installtest-defining-tests:

Defining Tests
--------------

Test cases are defined from TOML files.

The :repository-dir:`tests/install-tests/cases/manual_makefile` directory shows a fairly straight-forward example of a test.
Let's consider the configuration file from that directory:

.. literalinclude:: ../../tests/install-tests/cases/manual_makefile/conf.toml
   :language: toml

The ``[shared]`` table defines basic setup performed in all test cases enumerated by this file.
Each test-case is enumerated in a subtable of the ``test`` table (while ``[test]`` never appears in this file, it is implicitly defined because of the appearance of ``[test.<NAME>]``).
In this file, there are 2 tests, specified in the ``[test.simple]`` table and the ``test.SharedAndStatic`` table.
The testing program will refer to these tests as ``manual_makefile.simple`` and ``manual_makefile.SharedAndStatic``.

We go into detail about all available parameters down below:


.. embed-cli-output:: ../../tests/install-tests/installtest.py
   :args: entrydoc param


Design Rationale
----------------


Big Picture Ideas
+++++++++++++++++

While designing the test harness, we had the following (somewhat related) guiding principles:

- The tests should not require ``gracklepy`` to be installed.

  * Ideologically, these tests pertain to the core library and have nothing to do with ``gracklepy``

  * Practically, there's value to being able to run the tests without installing ``gracklepy``.
    Specifically, these tests explicitly check for the exact kinds of issues that would break builds of ``gracklepy``.

  * In the future, we may want to add an installtest checking the various approaches for installing ``gracklepy``.

- These tests should be portable and easy to run:

  * In practice, this means that the runner should use the system provided python installation.

  * currently, the tests just require the installation of docker (while that's a bit of a hurdle, it's :ref:`somewhat unavoidable <installtest-why-docker>`)

  * while I'm not opposed to using external python packages in the test harness, we would probably want to retain a separate pyproject.toml file (or at least a REQUIREMENTS.txt file) for the harness.

- Retain the capacity to easily integrate these tests with CTest

  * if the tests were "integrated," then CTest could be used execute individual test-cases or all of the tests.
    This could make the process of launching tests related to the core library more convenient.

  * for added clarity, all of our googletest unit tests are integrated with CTest (for that precise reason, the test harness's CLI takes cues from googletest's interface).

  * at present, we have not performed this integration, but we still want to retain the capacity to do it in the future.\ [#ctest-integration]_

Why don't we use ``pytest``?
++++++++++++++++++++++++++++

The short answer is: "it wouldn't really help us."

We provide several more detailed reasons:

1. ``pytest`` isn't really designed for these kinds of tests; it's primarily designed to test python code.

2. Since ``pytest`` is highly extensible, we could definitely extend it to perform these tests.
   In practice, we would have to write much of the code we already have (using ``pytest`` might realistically save us 100-200 lines of codes).

3. Using ``pytest`` comes with a few drawbacks:

   * Using ``pytest`` sacrifices the ability to easily integrate each installtest with CTest (this may not a big deal, but I'm a little hesitant to lock us into this choice, right now).

   * Refactoring the creation of docker images to use the fixture system would take some effort (it's definitely doable).

   * For non-experts, ``pytest``'s control flow is fairly tricky to follow.

   * Least importantly, adopting ``pytest`` would require us to determine how to track installtest.py's external dependencies (since ``pytest`` would be installtest.py's very first external dependnecy).


.. _installtest-why-toml:

Why we use TOML files for configuring tests
+++++++++++++++++++++++++++++++++++++++++++

To start, we decided not to define each test in pure python because it's convenient to define the tests in the same directory as the build-configuration files that are used in the tests.
While we could put a conf.py file in each directory, it's not obvious how to make python load multiple files named conf.py (but I'm confident that it's possible).

This leaves us with picking from a configuration file format.
The primary contenders are INI, JSON, TOML, and YAML.

**Why not INI?**
In short, INI is less well defined than any other format.
The validity of syntax in an INI depends on the parser, while a common standard defines valid syntax in all other cases.

**Why not JSON?**
We decided against JSON in order to be able to embed explanatory comments.
For context, the JSON specification does not support comments.

**Why TOML over YAML?**
TOML is preferable to YAML for a few reasons:

1. TOML is "just" sophisticated enough to meet our needs.
   While YAML is significantly more powerful, it is also significantly more complex and produces subtle edge cases.
   Some examples are listed `here <https://hitchdev.com/strictyaml/why/implicit-typing-removed/>`__ .

2. TOML is also preferable from a dependency perspective:

   - Currently, we include :source:`tests/install-tests/vendored_tomli.py` which is a vendored version of the `tomli <https://pypi.org/project/tomli/>`__ package (the file is ~800 lines).

  - Although vendoring is obviously undesirable, it's very much a short term solution.
    Once we are willing to assume everybody has a minimum python version of 3.11, we can just use python's builtin ``tomllib`` package.

  - In comparison, vendoring ``PyYAML`` is simply out of the question


.. _installtest-why-docker:

Why use Docker?
+++++++++++++++

We use docker for several reasons:

1. Most importantly, sandboxing the tests within a container is EXTREMELY convenient:

   * It lets tests make fairly invasive changes to the container, such as

     - installing software (e.g. ``pkg-config`` or a particular version of ``cmake``)

     - installing Grackle itself

     - installing a test-program

     - Overriding the ``LD_LIBRARY_PATH`` variable.

  * relatedly, it makes handling of dependencies a lot easier (we can simply install a C++ compiler or hdf5 ourselves rather than play games with a developer's setup).

2. It facillitates "caching" in 2 (related) ways.
   Since lots of tests are run after performing a common set of steps (e.g. installing to a local install-directory), we can build an image that stores the result of these steps, and then reuse these steps in multiple tests.
   Furthermore, we can also take advantage of built-in layer caching.

3. Finally, by using Docker, it becomes easy to run the tests on macOS (via Docker Desktop).

By using docker in these tests, we may also eventually see some benefits when it's time to start adding GPU support (we'll probably want to perform test builds in docker containers with CUDA/HIP when developing on platforms without supported GPUs).

Finally it's worth discussing a hypothetical concern: there may be drawbacks to tying these tests to a single external tool.
We aren't very worried because:

- If this ever does become an actual problem, we can always add support for using ``podman``, which was designed with the same interface as docker.

- In fact, the `cibuildwheel <https://cibuildwheel.pypa.io/en/stable/>`_ project demonstrated that it's pretty straight-forward to support both tools (the internal logic of installtest.py actually reuses some logic from cibuildwheel).



.. rubric:: Footnotes

.. [#packagers] It's unlikely that a packager would volunteer to maintain a Grackle package.
   This would require: (i) a package maintainers who themselves want to use Grackle (and want to contribute their time to the project), (ii) popular demand by lots of users of a package manager, or (iii) the use of Grackle as a dependency of other popular packages.
   None of these scenarios is likely.
   In any case, there may be some friction to doing this since Grackle isn't currently ABI-stable.

.. [#headeronly-lib] Header-only libraries don't necessarily need to be embedded in the source directory of downstream applications.
   But many of the discussed drawbacks are still relevant.

.. [#link-test] It *is* possible to temporarily overwrite the install-directory to try to temporarily install Grackle to that directory (by invoking ``cmake --install --prefix ``<DEST>``) and then test linking against the contents of that directory.
   However there are 2 wrinkles to that approach.
   First, not all install-tests will work very well in this case (think about install-tests that try to link a test-program against Grackle without using CMake or pkg-config).
   Second (and more importantly), different CMake logic is executed when the install destination is specified at initial configuration vs using the temporary overwrite logic (we have in fact had bugs pertaining to this logic).

.. [#ctest-integration] To integrate with CTest, we would need to make sure that we properly label each installtest test case and possibly disable them when a simple ``ctest`` command is invoked since these tests take a lot longer than the other kinds of tests (this is doable, but would take a little effort).
   We would also need to properly record each test case's dependency on the task of building a docker image (again, very doable).
