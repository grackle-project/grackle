.. _testing:

Testing
=======

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

.. _test_without_answer_verification:

Tests Without Answer Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you only want to quickly verify that everything runs, you can skip
generating the answer test results using the latest gold standard.

Once you have installed :ref:`pygrackle and the development dependencies <pygrackle-dev>`, this is relatively straight-forward.

Simply invoke ``py.test`` from the root of the Grackle directory with the ``--answer-dir=<PATH/TO/ANSWER-DIR>`` option **and** the ``--answer-store`` option.
The former option lets you can specify an arbitrary path where the test-answers will be recorded.
A concrete example of what this looks like is shown down below.

.. code-block:: shell-session

   ~ $ cd grackle
   ~/grackle $ py.test --answer-dir=./my_test_answers --answer-store

The above snippet instructs each answer-test to store the test result in the directory called **./my_test_answers**; an answer-test reports that it has "passed" as long as it is able to successfully store the result.
After you run the test-suite, you can delete the **./my_test_answers** directory (since we don't actually care about the test results).

When you launch the tests, the output will look like the following:

.. code-block:: shell-session

  ============================= test session starts ==============================
  platform darwin -- Python 3.11.9, pytest-8.2.1, pluggy-1.5.0
  rootdir: /Users/britton/Documents/work/research/simulation/grackle/grackle-git/src/python
  configfile: pyproject.toml
  testpaths: src/python/tests
  plugins: cov-5.0.0
  collected 62 items

  src/python/tests/test_chemistry.py ....                                  [  6%]
  src/python/tests/test_chemistry_struct_synched.py .                      [  8%]
  src/python/tests/test_dynamic_api.py ...                                 [ 12%]
  src/python/tests/test_get_grackle_version.py .                           [ 14%]
  src/python/tests/test_initialisation.py .                                [ 16%]
  src/python/tests/test_local_functions.py .                               [ 17%]
  src/python/tests/test_models.py .......................................  [ 80%]
  src/python/tests/test_primordial.py .                                    [ 82%]
  src/python/tests/test_query_units.py ...                                 [ 87%]
  src/python/tests/test_specific_heating_rate.py ....                      [ 93%]
  src/python/tests/test_volumetric_heating_rate.py ....                    [100%]

  ======================== 62 passed in 95.46s (0:01:35) =========================

Now it's time to :ref:`integrate grackle into your simulation code
<integration>`.

Tests With Answer Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To include answer test verification, one must first generate the
answers from the last gold standard version of Grackle. Gold standards
are marked with annotated git tags and are named 'gold-standard-v' and
then a number, for example, 'gold-standard-v1'. To find the latest
gold standard, type ``git tag``. If nothing is output, you may need to
first fetch the tags from the main repo with something like ``git
fetch origin --tags`` (where 'origin' in this example is assumed to
point to the main Grackle repository on github.)

To generate test results from the latest gold standard, follow these
steps:

#. Navigate your Grackle repository to the gold standard commit with,
   for example, ``git checkout gold-standard-v1``.

#. Re-compile the Grackle library and :ref:`re-install pygrackle
   <install-pygrackle>`.

#. Execute the test suite using command-line flags to instruct the test-runner to run answer tests in *store-mode*.
   This was already illustrated :ref:`above <test_without_answer_verification>`, we repeat the instructions here:

   - Execute ``py.test`` from the **src/python** directory while specifing the path  to the *answer-directory* with ``--answer-dir=<PATH/TO/ANSWER-DIR>`` **and** specifying the ``--answer-store`` option.

   - Be aware, this will directly overwrite any files that were previously stored in the answer-dir.

   - If we wanted to store the test-answers in **./my-test-answers**, we would invoke:

     .. code-block:: shell-session

        ~/grackle $ py.test --answer-dir=./my_test_answers --answer-store

#. Return to the branch of the repository you started with. If you just
   cloned the main repository, this will be called 'main', in which
   case do ``git checkout main``.

#. Re-compile the Grackle library and :ref:`re-install pygrackle
   <install-pygrackle>`.

#. Run the test suite again using a command-line flag to instruct the test runner to evaluate the answer-tests in *compare-mode*:

   - you just need to specify the *answer-directory* path with ``--answer-dir=<PATH/TO/ANSWER-DIR>``.

   - Do **NOT** pass the ``--answer-store`` flag.
     The absence of this flag is how the test-runner knows to use *compare-mode*.
     (If the flag were present, then the answer-tests would overwrite the answers)

   - To compare against the results previously written to **./my-test-answers**, you would invoke:

     .. code-block:: shell-session

        ~/grackle $ py.test --answer-dir=./my_test_answers
