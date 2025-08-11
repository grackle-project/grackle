
.. _manage-data-files:

Datafile Management
===================

We provide a command line tool to optionally manage Grackle's datafiles.

At a Quick Glance
-----------------

Currently, this command line tool is only accessible when :ref:`pygrackle is installed <install-pygrackle>`.
To execute the tool execute

.. code-block:: shell-session

   $ python -m pygrackle <args>...

Where ``<args>...`` is replaced with one or more command-line arguments.
For example, ``fetch`` will invoke a subcommand that downloads all associated files (if they aren't already downloaded).
You can use the ``--help`` option to get a list of all subcommands.
You can also pass the ``--help`` option after the name of a subcommand (e.g. you can use ``fetch --help``) to get more details about subcommand-specific options.

.. note::

   At the moment, this functionality is most useful for pygrackle.
   In the near future [#df1]_\ , it will be possible install pygrackle without manually downloading the grackle repository.
   At that time, this will be the most efficient way to retrieve the files.
   The pygrackle examples and some of the pygrackle tests rely upon this functionality.
   However, you are free to completely ignore this functionality for your own purposes.

   There is ongoing work to implement functionality for the Grackle C library to directly access the datafiles managed by this tool.
   When these efforts are finished, we plan to additionally provide this command-line-tool as a standalone program that is always installed alongside Grackle (so that you can access this functionality without installing pygrackle)

Description
-----------

.. include:: ../../src/python/pygrackle/utilities/grdata.py
   :start-after: [[[BEGIN-SECTION:DESCRIPTION]]]
   :end-before: [[[END-SECTION:DESCRIPTION]]]



Motivation
----------

.. include:: ../../src/python/pygrackle/utilities/grdata.py
   :start-after: [[[BEGIN-SECTION:MOTIVATION]]]
   :end-before: [[[END-SECTION:MOTIVATION]]]


How it works
------------

.. include:: ../../src/python/pygrackle/utilities/grdata.py
   :start-after: [[[BEGIN-SECTION:INTERNALS-OVERVIEW]]]
   :end-before: [[[END-SECTION:INTERNALS-OVERVIEW]]]


Sample Directory Structure
++++++++++++++++++++++++++

Down below, we sketch out what the directory-structure might look like:


.. literalinclude:: ../../src/python/pygrackle/utilities/grdata.py
   :language: none
   :start-after: [[[BEGIN:DIRECTORY-CARTOON]]]
   :end-before: [[[END:DIRECTORY-CARTOON]]]


.. rubric:: Footnotes

.. [#df1] Once `GH-#208 <https://github.com/grackle-project/grackle/pull/208>`__ is merged, you will be able to instruct pip to install pygrackle by just specifying the URL of the GitHub repository.
          We also have plans to upload pygrackle to pip.
