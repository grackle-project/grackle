
.. _examples:

Example Executables
===================

The Grackle repository provides example C, C++, and Fortran code for interacting with Grackle.

Example Descriptions
--------------------

The examples are located in the **src/example** directory.
The examples include:

    * **c_example.c** - C example

    * **c_local_example.c** - C example using only :ref:`local_functions`

    * **cxx_example.C** - C++ example

    * **cxx_omp_example.C** - C++ example using OpenMP

    * **fortran_example.F** - Fortran example

.. _how-to-run-example:

Preparing and Executing the Examples
------------------------------------

In this section, we explain how to prepare and execute the example executables.
Running an example is a useful way to quickly check whether Grackle is functioning correctly
(to more rigorously check that Grackle is fully functional, you can try :ref:`running the
test suite <testing>`).

The instructions for building and executing the examples vary based on the build-system.
In both cases, the examples require that you haven't cleanup up from your build.

1. Compile the example:

   .. tabs::

      .. group-tab:: Classic Build System
     
         Once you have already installed the grackle library, you can build the examples by typing ``make`` and the name of the file without extension.
         Assuming that you were just in the **src/clib** subdirectory, you would execute the following to build the C++ example:

         .. code-block:: shell-session

            ~/grackle/src/clib $ cd ../example
            ~/grackle/src/example $ make clean 
            ~/grackle/src/example $ make 

            Compiling cxx_example.C
            Linking
            Success!

      .. group-tab:: CMake Build System
 
         By default, the examples are automatically built with the rest of Grackle.
         The compiled example binaries can be found within **<build-dir>/example**, where **<build-dir>** is the arbitrary build-directory that you need to specify when compiling Grackle.

         .. warning::

            It's important that **<build-dir>** is a top-level directory in the grackle repository (e.g. something like **~/grackle/my-build** is fine, but choices like **~/grackle/../my-grackle-build** and **~/grackle/my_builds/my-first-build** are problematic).
            If this isn't the case, then the examples won't be able to locate the input data files.

   .. important::

      If you're using the Classic build system, make sure to add the path to the directory containing the installed **libgrackle.so** to your LD_LIBRARY_PATH (or DYLD_LIBRARY_PATH on Mac).
      This is **NOT** necessary for the CMake build system.
      More information is provided below.\ [#f1]_

2. Now we execute the example

   .. note::

      The examples make certain assumptions about the location of the data files.
      To ensure that the data files can be found, you should execute each example-binary from the same directory where the example binary is produced.

   To execute the example, invoke:

   .. tabs::

      .. group-tab:: Classic Build System

         .. code-block:: shell-session

            ~/grackle/src/example $ ./cxx_example

      .. group-tab:: CMake Build System

         .. code-block:: shell-session

            ~/grackle $ cd <build-dir>/examples
            ~/grackle/<build-dir>/examples $ ./cxx_example

   The output will look like the following:

   .. code-block:: shell-session

    The Grackle Version 2.2
    Mercurial Branch   default
    Mercurial Revision b4650914153d

    Initializing grackle data.
    with_radiative_cooling: 1.
    primordial_chemistry: 3.
    metal_cooling: 1.
    UVbackground: 1.
    Initializing Cloudy cooling: Metals.
    cloudy_table_file: ../../input/CloudyData_UVB=HM2012.h5.
    Cloudy cooling grid rank: 3.
    Cloudy cooling grid dimensions: 29 26 161.
    Parameter1: -10 to 4 (29 steps).
    Parameter2: 0 to 14.849 (26 steps).
    Temperature: 1 to 9 (161 steps).
    Reading Cloudy Cooling dataset.
    Reading Cloudy Heating dataset.
    Initializing UV background.
    Reading UV background data from ../../input/CloudyData_UVB=HM2012.h5.
    UV background information:
    Haardt & Madau (2012, ApJ, 746, 125) [Galaxies & Quasars]
    z_min =  0.000
    z_max = 15.130
    Setting UVbackground_redshift_on to 15.130000.
    Setting UVbackground_redshift_off to 0.000000.
    Cooling time = -1.434987e+13 s.
    Temperature = 4.637034e+02 K.
    Pressure = 3.345738e+34.
    gamma = 1.666645e+00.


.. rubric:: Footnotes

.. [#f1] In more detail, both build-systems use copies of the grackle-library within the build directory while compiling the example.

   * the Classic build-system **always** links Grackle against the shared-library version of Grackle and requires that Grackle is fully installed in a location known by the system (either a standard system location OR a location specified by ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH``).
   * In contrast, cmake automatically takes special-steps to try to ensure that each example-binary will link to the copy of the Grackle library (whether it is shared or static) that is in the ``<build-dir>``; in fact, Grackle doesn't even need to be installed to run the Grackle library.
   * With that said, if you compile Grackle as a shared library in a cmake build, an example-binary **might** try to use a copy of a shared grackle library found in a directory specified by ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH`` if one exists.
     The exact behavior may be platform dependent and also depends on whether CMake instructs the linker to use RPATH or RUNPATH (this is not specified by the cmake docs).

