.. _testing:

Running the Tests
=================

Grackle contains a number of unit and answer tests to verify that
everything is working properly.

The tests are primarily organized into 2 test suites:

1. the :ref:`core-library test suite <corelib-suite>`, which performs tests on the Core library
2. the :ref:`pygrackle test suite <pygrackle-suite>`, which performs tests on the pygrackle bindings

Our continuous integration system is also set up to ensure that all Python code conforms to `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`__

Historically, the pygrackle test suite included **all tests**.
More recently, the core-library test suite was introduced to make it easier to run unit tests on functionality that isn't directly exposed through pygrackle and to run tests involving compilation.
At this point, the pygrackle test suite includes a mix of unit tests and "answer tests" (we describe "answer tests" more down below). 

.. _corelib-suite:

Running the Core-Lib Test Suite
-------------------------------

The core-library test suite includes a mix of unit tests and "compilation tests."
The "compilation tests" verify that all code examples build, run, and return consistent results (these mostly check that the Grackle examples produce reasonable results and that the results are the same regardless of the choice of language/API bindings).

The core-library tests are driven by the ``ctest`` program that is shipped as part of ``CMake``.

- the ``CTest`` test-driver is integrated with ``CMake`` and it is designed to flexibly run various kinds of tests by invoking the command-line

- the unit tests that the ``CTest`` driver invokes are all implemented using the `GoogleTest testing Framework <https://google.github.io/googletest/>`__.
  Note: the dependency on GoogleTest is handled automatically (if you don't configure the build-system to use an existing installation, the build system will know to automatically fetch/compile the framework when instructed to compile the tests).

To configure Grackle builds to run the core-library tests, you must enable the ``GRACKLE_BUILD_TESTS`` option when you configure a CMake build of Grackle.
This option just enables ``CTest`` and some extra build-recipes needed to run the tests (it has no impact on how libgrackle is compiled/linked).

The following snippet illustrates how you might invoke the unit tests when performing a fresh build of grackle, with the :ref:`CMake build system <cmake_build>`, from the root of the Grackle repository:

.. code-block:: shell-session

   ~/grackle $ cmake -DGRACKLE_USE_DOUBLE=ON   \
                     -DGRACKLE_BUILD_TESTS=ON  \
                     -B<build-dir>
   ~/grackle $ cmake --build <build-dir>
   ~/grackle $ cd <build-dir>
   ~/grackle/<build-dir> $ ctest --output-on-failure

In the above snippet, ``<build-dir>`` should be replaced with your chosen build-directory's name (see the section on :ref:`CMake builds <cmake_build>` for more info); a common choice might be ``build``.
Let's quickly talk through the commands:

- First, we configure the build.
  You can add :ref:`other configuration options <available_cmake_options>` (e.g. ``GNinja``, ``GRACKLE_USE_OPENMP``, ``BUILD_SHARED_LIBS``, ``CMAKE_BUILD_TYPE``, etc.)
- Next, we invoke the build
- The final 2 commands move into the build-directory and invoke the tests from there (the ``--output-on-failure`` flag is completely optional).
  Starting in CMake 3.20 these 2 commands can be replaced with ``cmake --output-on-failure --test-dir <build-dir>``

When you launch the tests, the output will look like the following:

.. code-block:: shell-session

   Test project /Users/mabruzzo/packages/c++/grackle/build
         Start  1: GrExample.SatisfyPrereqs
    1/16 Test  #1: GrExample.SatisfyPrereqs .....................................   Passed    0.05 sec
         Start  2: GrExampleHistoricalCMP.c
    2/16 Test  #2: GrExampleHistoricalCMP.c .....................................   Passed    0.18 sec
         Start  3: GrExampleCMP.c_local
    3/16 Test  #3: GrExampleCMP.c_local .........................................   Passed    0.16 sec
         Start  4: GrExampleCMP.cxx
    4/16 Test  #4: GrExampleCMP.cxx .............................................   Passed    0.16 sec
         Start  5: GrExampleCMP.fortran
    5/16 Test  #5: GrExampleCMP.fortran .........................................   Passed    0.16 sec
         Start  6: InterpolationTest.Interpolate1D
    6/16 Test  #6: InterpolationTest.Interpolate1D ..............................   Passed    0.00 sec
         Start  7: InterpolationTest.Interpolate2D
    7/16 Test  #7: InterpolationTest.Interpolate2D ..............................   Passed    0.00 sec
         Start  8: InterpolationTest.Interpolate3D
    8/16 Test  #8: InterpolationTest.Interpolate3D ..............................   Passed    0.00 sec
         Start  9: InterpolationTest.Interpolate3Dz
    9/16 Test  #9: InterpolationTest.Interpolate3Dz .............................   Passed    0.00 sec
         Start 10: InterpolationTest.Interpolate2Df3D
   10/16 Test #10: InterpolationTest.Interpolate2Df3D ...........................   Passed    0.00 sec
         Start 11: InterpolationTest.Interpolate4D
   11/16 Test #11: InterpolationTest.Interpolate4D ..............................   Passed    0.00 sec
         Start 12: InterpolationTest.Interpolate5D
   12/16 Test #12: InterpolationTest.Interpolate5D ..............................   Passed    0.00 sec
         Start 13: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/0
   13/16 Test #13: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/0 ...   Passed    0.01 sec
         Start 14: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/1
   14/16 Test #14: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/1 ...   Passed    0.02 sec
         Start 15: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/2
   15/16 Test #15: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/2 ...   Passed    0.02 sec
         Start 16: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/3
   16/16 Test #16: VaryingPrimordialChem/APIConventionTest.GridZoneStartEnd/3 ...   Passed    0.02 sec

   100% tests passed, 0 tests failed out of 16

   Total Test time (real) =   0.80 sec

.. _pygrackle-suite:

Running the pygrackle Test Suite
--------------------------------

As already noted, the pygrackle suite includes unit tests and answer tests.

Unit tests (i.e., those with explicitly known correct answers) include
the following:

- correct library versioning

- correct behavior of the dynamic API

- proper and comoving unit systems are consistent

- mean molecular weight increases with metallicity

- atomic, primordial collisional ionization equilibrium agrees with
  the analytical solution

Answer tests are those whose correct answers must be generated from a
prior, trusted version of Grackle (i.e., the "gold standard"). The
tests are first run using this trusted version to generate the
results (in *store-mode*), then run again on the latest version to
compare (in *compare-mode*).
These tests include:

- all python examples run and give correct results for a range of
  parameter values

- all grackle 'calculate' functions return correct results for sets
  of random field values

We refer to the location where the results of answer-tests are stored as the "answer-directory." This is an arbitrary user-specified location.

.. note::

   The test of the **src/python/examples/yt_grackle.py** python-example requires that you supply test-data.
   The location of the test-data is specified through the :envvar:`!YT_DATA_DIR` environment variable.
   The **scripts/ci/fetch-test-data.py** script is provided as a convenience to fetch this data after the environment variable has been set.

Quick Primer on the Test Runner's CLI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the pytest test-runner always runs all available test cases.

- The suite's unit tests are **ALWAYS** available.

- By default, **all** answer-tests are fully disabled.
  These tests are made available, in *store-mode* or *compare-mode*, when the :option:`--answer-dir` command-line option is provided.
  The :option:`--answer-store` command line flag enables *store-mode*, while its absence enables *compare-mode*.

  .. option:: --answer-dir=<PATH>

     Specifies the path to the "answer-directory".
     This is the custom user-specified directory where answer-tests are stored (in *store-mode*) or read from (in *compare-mode*).

  .. option:: --answer-store

     The presence of this flag enables *store-mode*, where answer-tests results are stored to the answer-directory (any previously recorded answers will be overwritten).
     When :option:`--answer-dir` is specified and this flag is omitted, *compare-mode* is enabled.

*For contributors:* you may find pytest's `build-in command-line interface <https://docs.pytest.org/en/stable/how-to/usage.html>`__ useful during debugging (e.g. you can instruct pytest to only run a subset of all available tests).


Sample Usage of the pygrackle Test Suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following snippets illustrate ways how you might invoke the pygrackle test suite in different scenarios.

Unless noted otherwise, each scenario assumes that you have already :ref:`installed pygrackle <install-pygrackle>` (reminder, you currently need editable installs for these tests).
Each shows snippets assuming that you are at the root of the Grackle directory.

.. tabs::

   .. tab:: No Answer Tests

      In this scenario, we illustrate how to run the unit-tests in the pygrackle test suite and skip all answer tests (**NOTE:** there currently isn't an easy way to do the opposite).
      To do this, invoke:

      .. code-block:: shell-session

         ~/grackle $ py.test

      The output will look resemble the following:

      .. code-block:: shell-session

         ============================= test session starts ==============================
         platform darwin -- Python 3.10.14, pytest-8.3.5, pluggy-1.5.0
         rootdir: /Users/mabruzzo/packages/c++/grackle
         configfile: pyproject.toml
         testpaths: src/python/tests
         plugins: anyio-4.2.0
         collected 62 items

         src/python/tests/test_chemistry.py ....                                  [  6%]
         src/python/tests/test_chemistry_struct_synched.py .                      [  8%]
         src/python/tests/test_dynamic_api.py ...                                 [ 12%]
         src/python/tests/test_get_grackle_version.py .                           [ 14%]
         src/python/tests/test_initialisation.py s                                [ 16%]
         src/python/tests/test_local_functions.py s                               [ 17%]
         src/python/tests/test_models.py sssssssssssssssssssssssssssssssssssssss  [ 80%]
         src/python/tests/test_primordial.py .                                    [ 82%]
         src/python/tests/test_query_units.py ...                                 [ 87%]
         src/python/tests/test_specific_heating_rate.py ....                      [ 93%]
         src/python/tests/test_volumetric_heating_rate.py ....                    [100%]

         =========================== short test summary info ============================
         SKIPPED [1] src/python/tests/test_initialisation.py:103: no --answer-dir option found
         SKIPPED [1] src/python/tests/test_local_functions.py:191: no --answer-dir option found
         SKIPPED [39] src/python/tests/test_models.py:23: no --answer-dir option found
         ================== 21 passed, 41 skipped, 1 warning in 2.61s =================== 


   .. tab:: Full Suite (No Answer Verification)

      In this scenario, we illustrate how to quickly verify that all of the answer-test logic runs (without verifying that the test results are correct).

      For concreteness, we arbitrarily adopt an answer-directory named **./my_test_answers**, which we will delete at the end of the scenario (you probably don't want to use that name if the directory already exists).

      The following command runs the full test suite, with the answer-tests configured in *store-mode* (results are stored to **./my_test_answers**):

      .. include:: testing_snippets/generate.rst

      .. tip::

         After running the above, you can delete the answer-test directory (**./my_test_answers** in this scenario) since we don't actually care about the test results.


   .. tab:: Record Answers & then Run Full Suite

      This scenario comes up the very first time you want to run the answer tests (with verification) on a given machine.
      If you have already generated answers before, skip to the next tab.

      .. note::

         This scenario involves installing (and uninstalling) different versions of Pygrackle.
         These examples all use the simplest (and recommended) approach for installing Pygrackle, which builds it as a standalone library (and automatically manages all details about the core library).

         You could also run the tests using :ref:`the other methods of installing Pygrackle <install-pygrackle>`, but that involves extra steps (in those cases you need to manually build and install/link the associated version of the core library).
         

      The basic premise is that you want to run the full pygrackle test suite (with answer verification) on a target version of the code called ``<target>``.
      ``<target>`` is a placeholder for commit hash or a branch name (e.g. ``main``).

      .. rubric:: Basic Setup: Idetifying the Gold Standard Tag

      To include answer test verification, one must first generate the
      answers from the last gold standard version of Grackle. Gold standards
      are marked with annotated git tags and are named 'gold-standard-v' and
      then a number, for example, 'gold-standard-v1'. To find the latest
      gold standard, type ``git tag``. If nothing is output, you may need to
      first fetch the tags from the main repo with something like ``git
      fetch origin --tags`` (where 'origin' in this example is assumed to
      point to the main Grackle repository on GitHub.)

      .. rubric:: Generate Answers from Gold Standard

      For concreteness, let's assume that we are generating test results from ``gold-standard-v3`` and storing them in an answer-directory called **./my_test_answers**.

      We need to checkout the gold standard version and installs the included version of pygrackle:

      .. code-block:: shell-session

         ~/grackle $ git checkout gold-standard-v3 # replace <version>
         ~/grackle $ pip uninstall pygrackle
         ~/grackle $ pip install -e .

      The following command runs the full test suite, with the answer-tests configured in *store-mode* (results are stored to **./my_test_answers**):

      .. include:: testing_snippets/generate.rst

      .. rubric:: Run the full test suite on ``<target>``

      Now, it's time to go back to the ``<target>`` version. 
      The following checks out that version of Grackle and installs the included version of pygrackle:

      .. code-block:: shell-session

         ~/grackle $ git checkout <target>
         ~/grackle $ pip uninstall pygrackle
         ~/grackle $ pip install -e .

      Finally, the following command runs the full suite, with the answer-tests enabled in **compare-mode** (answer tests only pass if the results are consistent with the answers we previously recorded to the answer directory).

      .. code-block:: shell-session

         ~/grackle $ ~/grackle $ py.test --answer-dir=./my_test_answers

      The output of this command is in the same style we showed earlier

      .. warning::

         Do **NOT** pass the ``--answer-store`` flag.
         The absence of this flag is how the test-runner knows to use *compare-mode*.
         (If the flag were present, then the answer-tests would overwrite the answers)

   .. tab:: Run Full Suite (pre-recorded answers)

      This final scenario assumes that you previously recorded the answers to the tests from the gold-standard.
      For concreteness, let's assume that the answers were recorded to an answer-directory called **./my_test_answers**.

      The following command runs the full suite, with the answer-tests enabled in **compare-mode** (answer tests only pass if the results are consistent with the answers we previously recorded to the answer directory).

      .. code-block:: shell-session

         ~/grackle $ ~/grackle $ py.test --answer-dir=./my_test_answers

      

