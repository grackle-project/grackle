.. _testing:

Running the Test Suite
----------------------

Grackle contains a number of unit and answer tests to verify that
everything is working properly.

Unit tests (i.e., those with explicitly known correct answers) include
the following:

 - correct library versioning

 - correct behavior of the dynamic API

 - proper and comoving unit systems are consistent

 - mean molecular weight increases with metallicity

 - atomic, primordial collisional ionization equilibrium agrees with
   the analytical solution

 - all Python code conforms to `PEP 8
   <https://www.python.org/dev/peps/pep-0008/>`__

Answer tests are those whose correct answers must be generated from a
prior, trusted version of Grackle (i.e., the "gold standard"). The
tests are first run using this trusted version to generate the
results, then run again on the latest version to compate. These tests
include:

 - all code examples build, run, and return correct results

 - all python examples run and give correct results for a range of
   parameter values

 - all grackle 'calculate' functions return correct results for sets
   of random field values

.. _test_without_answer_verification:

Tests Without Answer Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you only want to quickly verify that everything runs, you can skip
generating the answer test results using the latest gold standard.

Once you have installed :ref:`pygrackle and the development dependencies <pygrackle-dev>`,
there are currently 2 ways to do this:

.. tabs::

   .. group-tab:: Command Line Flag

      Simply invoke ``py.test`` from the **src/python** directory with the ``--answer-store`` flag:

      .. code-block:: shell-session

         ~ $ cd grackle/src/python
         ~/grackle/src $ py.test --answer-store

      This will both invoke the test suite and tell it to simply generate new answers and not compare with previously existing ones.

   .. group-tab:: Environment Variable

      First set the ``GENERATE_PYGRACKLE_TEST_RESULTS`` environment variable to 1.

      .. code-block:: shell-session

          ~ $ export GENERATE_PYGRACKLE_TEST_RESULTS=1

      This will tell the test suite to simply generate new answers and not
      compare with previously existing ones.

      The tests can then be run from the **src** directory by typing ``make test``:

      .. code-block:: shell-session

         ~ $ cd grackle/src
         ~/grackle/src $ make test

      or from the **src/python** directory by typing ``py.test``:

      .. code-block:: shell-session

         ~ $ cd grackle/src/python
         ~/grackle/src $ py.test

      .. warning::

         The ability to use an environment variable may be removed.

Once you launch the test, the output will look like the following:

.. code-block:: shell-session

  ==================================== test session starts ====================================
  platform darwin -- Python 3.11.9, pytest-8.2.1, pluggy-1.5.0
  rootdir: /Users/britton/Documents/work/research/simulation/grackle/grackle-git/src
  plugins: cov-5.0.0
  collected 65 items

  python/tests/test_chemistry.py ....                                                   [  6%]
  python/tests/test_chemistry_struct_synched.py .                                       [  7%]
  python/tests/test_code_examples.py ......                                             [ 16%]
  python/tests/test_dynamic_api.py ...                                                  [ 21%]
  python/tests/test_get_grackle_version.py .                                            [ 23%]
  python/tests/test_initialisation.py .                                                 [ 24%]
  python/tests/test_local_functions.py .                                                [ 26%]
  python/tests/test_models.py .......................................                   [ 86%]
  python/tests/test_primordial.py .                                                     [ 87%]
  python/tests/test_specific_heating_rate.py ....                                       [ 93%]
  python/tests/test_volumetric_heating_rate.py ....                                     [100%]

  ============================== 65 passed in 148.47s (0:02:28) ===============================

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

#. Execute the test suite with instructions to generate test results.
   The 2 ways to do this are described :ref:`above <test_without_answer_verification>`: (i) execute ``py.test`` with the ``--answer-store`` flag or (ii) set the ``GENERATE_PYGRACKLE_TEST_RESULTS`` environment variable to 1 before executing the test suite.
   By default, this will create test result files in the directory **src/python/tests/test_answers**.

#. Return to the branch of the repository you started with. If you just
   cloned the main repository, this will be called 'main', in which
   case do ``git checkout main``.

#. Re-compile the Grackle library and :ref:`re-install pygrackle
   <install-pygrackle>`.

#. If you previously assigned a value to the ``GENERATE_PYGRACKLE_TEST_RESULTS`` variable, you must now unset the variable **OR** assign it a value of 0.

#. Run the test suite again (do **NOT** pass the ``--answer-store`` flags.
   This time, the answer tests will be compared with the previously generated results.

Other Test Configuration
^^^^^^^^^^^^^^^^^^^^^^^^
To run the test-suite without any answer tests, invoke ``py.test`` with the ``--answer-skip`` flag.

To control the location of the directory where the test answers are save, you can invoke ``py.test`` with the ``--local-dir`` flag


