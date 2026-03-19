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

All of these steps are completed within a ``docker container``.

How the tests work
------------------

.. todo::
   ADD ME

Design Decisions
----------------

.. todo::
   ADD ME

Running the tests
-----------------

.. todo::
   ADD ME


Adding a new test case
----------------------

.. todo::
   ADD ME

