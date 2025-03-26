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
results (in *store-mode*), then run again on the latest version to
compare (in *compare-mode*).
These tests include:

 - all code examples build, run, and return correct results

 - all python examples run and give correct results for a range of
   parameter values

 - all grackle 'calculate' functions return correct results for sets
   of random field values

We refer to the location where the results of answer-tests are stored as the "answer-directory." This is an arbitrary user-specified location.

Quick Primer on the Test Runner's CLI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, the pytest test-runner always runs all available test cases.

 - The unit tests are **ALWAYS** available.

 - By default, **all** answer-tests are fully disabled.
   Command-line options make them available in *store-mode* or *compare-mode*.
   The ``--answer-dir=<PATH/TO/ANSWER-DIR>`` flag is required for both modes; it specifies the path to the user's chosen "answer-directory" (where answer-tests results are stored/read-from).
   The ``--answer-store`` flag enables *store-mode*, while its absence enables *compare-mode*.

*For contributors:* you may find pytest's `build-in command-line interface <https://docs.pytest.org/en/stable/how-to/usage.html>`__ useful during debugging (e.g. you can instruct pytest to only run a subset of all available tests).

.. _test_without_answer_verification:

Tests Without Answer Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you only want to quickly verify that everything runs, you can skip
generating the answer test results using the latest gold standard.

Once you have installed :ref:`pygrackle and the development dependencies <pygrackle-dev>`, this is relatively straight-forward.

Simply invoke ``py.test`` from the **src/python** directory with the ``--answer-dir=<PATH/TO/ANSWER-DIR>`` option **and** the ``--answer-store`` option.
For the former You can specify an arbitrary path where the test-answers will be recorded.
A concrete example of what this looks like is shown down below.

.. code-block:: shell-session

   ~ $ cd grackle/src/python
   ~/grackle/src/python $ py.test --answer-dir=./my_test_answers --answer-store

The above snippet instructs each answer-test to store the test result in the directory called **./my_test_answers**; an answer-test reports that it has "passed" as long as it is able to successfully store the result.
After you run the test-suite, you can delete the **./my_test_answers** directory (since we don't actually care about the test results).

When you launch the tests, the output will look like the following:

.. code-block:: shell-session

  ==================================== test session starts ====================================
  platform darwin -- Python 3.11.9, pytest-8.2.1, pluggy-1.5.0
  rootdir: /Users/britton/Documents/work/research/simulation/grackle/grackle-git/src/python
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

#. Execute the test suite using command-line flags to instruct the test-runner to run answer tests in *store-mode*.
   This was already illustrated :ref:`above <test_without_answer_verification>`, we repeat the instructions here:

   - Execute ``py.test`` from the **src/python** directory while specifing the path  to the *answer-directory* with ``--answer-dir=<PATH/TO/ANSWER-DIR>`` **and** specifying the ``--answer-store`` option.

   - Be aware, this will directly overwrite any files that were previously stored in the answer-dir.

   - If we wanted to store the test-answers in **./my-test-answers**, we would invoke:

     .. code-block:: shell-session

        ~/grackle/src/python $ py.test --answer-dir=./my_test_answers --answer-store

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

        ~/grackle/src/python $ py.test --answer-dir=./my_test_answers
