.. _install-tests:

installtest.py
==============

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


Why don't other open-source projects place emphasis on these tests?
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


The reader may understandably wonder "why don't other open-source software library projects have such a big emphasis on these sorts of tests?"
To start, Grackle is a library, not an executable (the sorts of things we're testing generally aren't an issue for executables).

But, let's compare Grackle against other open-source libraries.
There are several reasons for our emphasis:

1. Developers of other open-source software libraries often place much less focus on making their projects easy-to-install.

   .. COMMENT
      We touch on this point (i.e. why we care more about describing the installation process and making it as easy as possible) more than other projects in :ref:`cmake_buildsystem_design_rationale`.
      If we feel compelled to touch on this point again, I really think we should relocate this explanation, and just link to it.

   * They commonly assume that the only people that will be building a software library are developers who are already well-versed in building/linking/consuming software libraries.

   * When popular libraries (think of hdf5, fftw3, libpng, zlib, xzutils) are used in projects written less experienced developers, the projects are typically made available by package managers (e.g. apt/dpkg, yum/rpm, brew).
     Essentially, package maintainers (rather than the project developers) take responsibility for distributing the program.
     At this time, no packager has expressed any interest in doing this for us.\ [#packagers]_

   * When a project does want to make a project easy to distribute, they common way to try to circumvent issues by distributing libraries as a header-only library or as a single source and header library (the idea is that downstream projects directly embed the dependency in their source code code repository).
     Setting aside the issues this might cause for getting improvements/bug fixes contributed back to the main repository, this isn't very viable for us.
     First, it only works for downstream projects written in the same standard of C++ as Grackle or newer (at present, C++ 17 or newer).
     Second, it won't really work be viable as we start introducing GPU support.
     The need for data files introduces additional complexity.

   * Another alternative is for a project to try to directly distribute pre-compiled copies of the library.
     This isn't a viable strategy since Grackle depends on HDF5, which historically hasn't been ABI-stable between versions.
     We can revisit this in the future, since starting with 2.0, it sounds like HDF5 developers have started ABI stability very seriously -- we need to wait for these releases to be widely used.\ [#h5-abi]_
     With that said, this isn't won't really work very well when we start to add support for different GPU backends since the GPU runtime libraries aren't necessarily ABI stable (plus we would need to distribute a separate version for each backend).

   In contrast, we assume that users of downstream codes probably need to manually install Grackle themselves, and they lack a lot of experience pertaining to the process of compiling/linking software.
   In fact, Grackle could be the very piece of code an undergrad physics/astronomy major or first-year grad student could compile.
   Thus, we make an effort to provide lots of instruction.


2. We attempted to make Grackle easy to use as a static library in pure C programs.
   It turns out that this is at the root of a lot of complexity:

   * complications arise because Grackle itself has external dependencies

     - it explicitly depends on hdf5.

     - it also has implicit dependencies a regular C compiler (or linker) won't know about out of the box:

       * Most obviously, this is the C++ standard library (historically, it was the Fortran standard library).

       * When compiled with OpenMP, Grackle has an additional implicit dependency on the library.

   * The fact that we also took steps to make it easy to install both shared and static copies of Grackle to a single directly adds additional complexity.

   * Moreover, we added slightly more complexity by making it possible for external CMake projects to more easily select between the use of shared/static libraries.

   .. note::

      To be clear, the only real alternative alternative to this choice of supporting Grackle as a static libraries is tell people that they are on their own outside of embedded Grackle builds (within external CMake projects).

3. Grackle doesn't have a stable ABI and we attempt to remain somewhat backwards compatible with choices made by the classic build-system.
   In practice, this translates to creating a ``libgrackle.so`` (or ``libgrackle.dylib`` file on MacOS) symlink during installation that links to a somewhat atypically named shared library during the install step (the source shared library is named to communicate libgrackle ABI incompatability).
   While the logic for doing this is simple, it needs testing (it has held bugs before).


Why not run these tests as part of the core test suite?
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

At a surface level, the steps of a single test-case makes it sound like these tests could be part of the core-library test suite.
In fact, the core-lib tests already involve a set of "compile-tests" that closely resemble the structure of some tests.
Moreover, other open source projects might try to test this sort of thing by constructing a ctest test-case that explicitly calls the ``ctest`` executable with the ``--build-and-test`` option (an example can be found `here <https://github.com/fmtlib/fmt/blob/dc05bee30755bb993add7ed90845d87c2315f9c6/test/CMakeLists.txt#L166>`__).

There are 2 limitation to this approach:

1. it **only** tests the direct product of a build (if you build Grackle as a shared library, you only test builds against a shared library)

2. CMake/CTest is designed for tests to occur when targets are found in the build-directly.
   In practice, most people link against Grackle after it has been fully installed.\ [#link-test]_

These limitations aren't usually as much of a concern in external project because (as noted above) there is usually far less emphasis on making the project easy to install/link against, for less experienced developers.


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
Each docker image 






When you invoke

In practice, step 1 consists of creating a docker image, while steps 2-4 are executed in a docker container based on the image.
The use of docker is intended to take adva

.. _installtest-docker-images:

Docker Images
-------------

As we've noted in multiple other places, each individual install-test is run inside of a Docker Container.
We provide further justification for this choice :ref:`here <

.. _installtest-defining-tests:

Defining Tests
--------------

.. embed-cli-output:: ../../tests/install-tests/installtest.py
   :args: entrydoc param


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

  * at present, we have not performed this integration, but we still want to retain the capacity to do it in the future.\ [#ctest-integration]_

Why don't we use pytest?
++++++++++++++++++++++++

The short answer is: "it wouldn't really help us."

We provide several more detailed reasons:

1. pytest isn't really designed for these kinds of tests; it's primarily designed to test python code.

2. Since pytest is extendible, we could definitely extend it to make it work for this purpose.
   In practice, we would have to write much of the code we already have (using pytest might realistically save us 100-200 lines of codes).

3. Using pytest comes with a few drawbacks:

   * Using pytest sacrifices the ability to easily integrate the install-tests with CTest (this may not a big deal, but I'm a little hesitant to lock us into this choice, right now).

   * Refactoring the creation of docker images to use the fixture system would take some effort (it's definitely doable).

   * For non-experts, pytest's control flow is fairly tricky to follow.

   * Least importantly, as the installtest.py's first external python dependency, using it would require us to give some thought about specifying python dependencies (not a big deal)


.. _installtest-why-toml:

Why we use TOML files for configuring tests
+++++++++++++++++++++++++++++++++++++++++++

To start, we decided not to define each test in pure python because it's convenient to define the tests in the same directory as the build-configuration files that are used in the tests.
While we could put a conf.py file in each directory, importing the logic gets a little tricky...

This leaves us with picking from a configuration file format.
The primary contenders are INI, JSON, TOML, and YAML.

**Why not INI?**
In short, INI is less well defined than any other format.
The validity of syntax in an INI depends on the parser, while a common standard defines valid syntax in all other cases.

**Why not JSON?**
We decided against JSON in order to be able to embed explanatory comments.

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
We aren't very worried because

- If this ever does become an actual problem, we can always add support for using ``podman``, which was designed with the same interface as docker.
- In fact, the `cibuildwheel <https://cibuildwheel.pypa.io/en/stable/>`_ project demonstrated that it's pretty straight-forward to support both tools (the internal logic of installtest.py actually reuses some logic from cibuildwheel).



.. rubric:: Footnotes

.. [#packagers] It's unlikely that a packager would volunteer to maintain a Grackle package.
   This would require: (i) a package maintainers who themselves want to use Grackle (and want to contribute their time to the project), (ii) popular demand by lots of users of a package manager, or (iii) the use of Grackle as a dependency of other popular packages.
   None of these scenarios is likely.
   In any case, there may be some friction to doing this since Grackle isn't currently ABI-stable.

.. [#h5-abi] If we reach a point where we are confident that most people are using versions of hdf5 that are 2.0 (or are compatible with version 2.0), we could compile and distribute a shared-lib version of Grackle that was linked against version 2.0 (i.e. the earliest version) and it *should* work on any downstream system (barring naming errors).

.. [#link-test] It *is* possible to temporarily overwrite the install-directory to try to temporarily install Grackle to that directory (by invoking ``cmake --install --prefix ``<DEST>``) and then test linking against the contents of that directory.
   However there are 2 wrinkles to that approach.
   First, not all install-tests will work very well in this case (think about install-tests that try to link a test-program against Grackle without using CMake or pkg-config).
   Second (and more importantly), different CMake logic is executed when the install destination is specified at initial configuration vs using the temporary overwrite logic (we have in fact had bugs pertaining to this logic).

.. [#ctest-integration] To integrate with CTest, we would need to make sure that we properly label each installtest test case and possibly disable them when a simple ``ctest`` command is invoked since these tests take a lot longer than the other kinds of tests (this is doable, but would take a little effort).
   We would also need to properly record each test case's dependency on the task of building a docker image (again, very doable).
