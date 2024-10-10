
.. _manage-data-files:

Datafile Management
===================

We provide a command line tool to optionally manage Grackle's datafiles.
We call this the ``grdata`` tool.

At a Quick Glance
-----------------

We provide 2 ways to access this tool:

1. As a standalone command line application installed alongside of grackle.
2. As a command line tool shipped :ref:`as a part of pygrackle <install-pygrackle>`.

To execute the tool:

.. tabs::

   .. group-tab:: As a Standalone CLI

      When you build grackle as a standalone application with the CMake build-system, the program can generally be found inside the build directory at *<path-to-build-dir>/grackle/bin/grdata*. [#df1]_
      Thus, you can invoke the tool with:

      .. code-block:: shell-session

         $ ./<path-to-build-dir>/grackle/bin/grdata <args>...

      When you install grackle (with either the traditional build-system OR the CMake build-system), the program will be installed alongside the core Grackle library.
      If the core Grackle library is at *<install-prefix>/lib/libgrackle...*, then the tool is probably installed to *<install-prefix>/bin/grdata*. [#df2]_
      In this case, you can invoke the tool with:


      .. code-block:: shell-session

         $ <install-prefix>/bin/grdata <args>...


      .. note::

         In detail, the tool is implemented as an executable python script (it only uses the python standard library) and it relies upon a `shebang <https://en.wikipedia.org/wiki/Shebang_%28Unix%29>`__ to launch the program.

         In the event that, that your machine finds an invalid python version (at this time of writing, the minimum required version is 3.6.1) or can't find any python interpretter, you have 2 options:

         1. you can modify the shebang path OR

         2. you can directly invoke the tool with python: ``<path/to/python> <path/to/grdata> <args>...``.


   .. group-tab:: As a part of Pygrackle

      When Pygrackle is installed, you can invoke

      .. code-block:: shell-session

         $ python -m pygrackle <args>...

      .. note::

         If you choose to install pygrackle without manually downloading the grackle repository, this tool is the most efficient way to download the files.



In the sample snippets ``<args>...`` is replaced with one or more command-line arguments.
For example, ``fetch`` will invoke a subcommand that downloads all associated files (if they aren't already downloaded).
You can use the ``--help`` option to get a list of all subcommands.
You can also pass the ``--help`` option after the name of a subcommand (e.g. you can use ``fetch --help``) to get more details about subcommand-specific options.


The pygrackle examples and the pygrackle tests all rely upon this functionality.
The Grackle C library has support for access the datafiles managed by this tool.
Some of the examples may soon rely upon the functionality.

.. important::

   Instances of the grdata tool are associated with a single version of Grackle (if you are using Pygrackle, the version of the core Grackle c-library is the relevant version number).

   Internal (Py)Grackle logic that automatically locate and use the files managed by this tool, will **only** work if the files have been "fetched" by versions of the grdata tool that **exactly** match.
   In detail, matching version strings have identical major, minor, and micro version numbers.
   Both version strings must either have an identical suffix (namely "-dev") or they must both lack a suffix.
   For development versions of Grackle, we do **NOT** require the commit-hashes to exactly match.

   Consider the following hypothetical scenario:

   * You install version 3.4.0 of Grackle and use the associated version of ``grdata`` to install datafiles.

      * The version of Pygrackle that wraps 3.4.0 will also be able to automatically locate these datafiles

      * If you install a different version of Grackle (say version 3.3.1-dev or 3.4.0-dev or 3.4.1 or 3.4.1-dev) or a copy of Pygrackle that wraps a different version, you will need to use the copy of ``grdata`` associated with that copy to fetch the datafiles.

   * You install the version of Pygrackle that wraps version 3.4.0 of Grackle and use it to install the associated datafiles

      * Any build/installation of Grackle version 3.4.0 (whether or not it is wrapped by Pygrackle) will also be able to locate and use these datafiles 
 
   As we will discuss, the grdata tool takes a steps to deduplicate the storage used by the datafiles for different Grackle versions.
   For example, if you use ``grdata`` instances to associate 5 distinct grackle versions only one copy of the **CloudyData_UVB=HM2012_high_density.h5** file's data will be stored on disk (i.e. that data will only take ~6.74 MB of disk space rather than ~33.7 MB)

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

.. [#df1] The precise location is sensitive to CMake's standard `CMAKE_INSTALL_BINDIR` variable.
   On just about all platforms this defaults to *bin*, but it is plausible that it could default to a different value.
   You are also free to override the value of the `CMAKE_INSTALL_BINDIR` variable.

.. [#df2] Again, the precise location of ``grdata`` installed by the CMake build-system is sensitive to CMake's standard `CMAKE_INSTALL_BINDIR` variable.
   Additionally the precise location of the core Grackle library is sensitive to the CMake's standard `CMAKE_INSTALL_LIBDIR` variable, which defaults to platform's preferred default values (again, you can overwrite it).
   CMake's `CMAKE_INSTALL_LIBDIR` variable is known to vary considerably more between platforms.
   Common defaults include *lib* (most common), *lib64*, *lib32*, or *libx32*.

