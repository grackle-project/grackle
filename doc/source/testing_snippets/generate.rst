.. code-block:: shell-session

   ~/grackle $ py.test --answer-dir=./my_test_answers --answer-store

Each answer-test executed by this command reports that it has "passed" as long as it is able successfully store the result.
The output will look like the following:

.. code-block:: shell-session

  ============================= test session starts ==============================
  platform darwin -- Python 3.11.9, pytest-8.2.1, pluggy-1.5.0
  rootdir: /Users/britton/Documents/work/research/simulation/grackle/grackle-git/python
  configfile: pyproject.toml
  testpaths: python/tests
  plugins: cov-5.0.0
  collected 62 items

  python/tests/test_chemistry.py ....                                  [  6%]
  python/tests/test_chemistry_struct_synched.py .                      [  8%]
  python/tests/test_dynamic_api.py ...                                 [ 12%]
  python/tests/test_get_grackle_version.py .                           [ 14%]
  python/tests/test_initialisation.py .                                [ 16%]
  python/tests/test_local_functions.py .                               [ 17%]
  python/tests/test_models.py .......................................  [ 80%]
  python/tests/test_primordial.py .                                    [ 82%]
  python/tests/test_query_units.py ...                                 [ 87%]
  python/tests/test_specific_heating_rate.py ....                      [ 93%]
  python/tests/test_volumetric_heating_rate.py ....                    [100%]

  ======================== 62 passed in 95.46s (0:01:35) =========================

