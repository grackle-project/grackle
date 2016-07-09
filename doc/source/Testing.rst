.. _testing:

Running the Test Suite
----------------------

Grackle contains a number of unit and answer tests to verify that everything is
working properly.  These will verify that:

 - proper and comoving unit systems are consistent

 - atomic, primordial collisional ionization equilibrium agrees with
   the analytical solution

 - all code examples build and run

 - all python examples run and give correct results

 - all Python code conforms to `PEP 8
   <https://www.python.org/dev/peps/pep-0008/>`__

Once you have installed :ref:`pygrackle <python>`, the tests can be run from the
**src** directory by typing ``make test``:

.. code-block:: bash

  ~ $ cd grackle/src
  ~/grackle/src $ make test

or from the **src/python** directory by typing ``py.test``:

.. code-block:: bash

  ~ $ cd grackle/src/python
  ~/grackle/src $ py.test

  ===================================== test session starts ======================================
  platform darwin -- Python 2.7.11, pytest-2.8.1, py-1.4.30, pluggy-0.3.1
  rootdir: /Users/britton/Documents/work/simulation/grackle/grackle/src/python, inifile:
  collected 13 items

  tests/test_chemistry.py ...
  tests/test_code_examples.py ....
  tests/test_examples.py ........
  tests/test_flake8.py .
  tests/test_primordial.py .

  ================================== 17 passed in 68.83 seconds ==================================

Now it's time to :ref:`integrate grackle into your simulation code
<integration>`.