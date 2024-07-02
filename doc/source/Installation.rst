.. _obtaining_and_building_enzo:

Installation
============

There are 3 steps to setting up Grackle on your system

   1. :ref:`Install Grackle's Dependencies <install_grackle_dependencies>`

   2. :ref:`Download Grackle <download_grackle>`

   3. Build and install Grackle using the :ref:`classic build system <classic_build>` or the :ref:`CMake build system <cmake_build>`.


.. _install_grackle_dependencies:

Dependencies
------------

In addition to C/C++ and Fortran compilers, the following dependency must 
also be installed:

   * `HDF5 <http://www.hdfgroup.org/HDF5/>`_, the hierarchical data format.
     HDF5 also may require the szip and zlib libraries, which can be
     found at the HDF5 website.  Compiling with HDF5 1.8 or greater
     requires that the compiler directive ``H5_USE_16_API`` be specified.
     This can be done with ``-DH5_USE_16_API``, which is in the machine 
     specific make files.

Although many systems already have them installed, both build systems have additional dependencies:

   * the :ref:`classic build system <classic_build>`, employs the ``makedepend`` and the `libtool <https://www.gnu.org/software/libtool/>`_ utilities.
     It's often easiest to download these dependencies through your system's package manager.

   * the :ref:`CMake build system <cmake_build>` requires cmake to be installed.
     It's easiest download a binary distribution from the `CMake website <https://cmake.org/download/>`_ or use your system's package manager.

.. _download_grackle:

Downloading
-----------

Grackle is available in a git repository
`here <https://github.com/grackle-project/grackle>`__. Excellent guides
to git and GitHub are available at
`guides.github.com <https://guides.github.com/>`__. To clone the Grackle
repo, do the following:

.. highlight:: none

::

    ~ $ git clone https://github.com/grackle-project/grackle

Additional files containing cooling tables and test results are stored in
a submodule linked to the Grackle repository. To get these, run the
following command from anywhere within the repository:

.. highlight:: none

::

    ~ $ git submodule update --init


.. _classic_build:

Building with Classic Build-System
----------------------------------

The classic compilation process for grackle is very similar to that of
`Enzo <http://enzo-project.org>`_.  For more details on the Enzo build 
system, see the `Enzo build documentation 
<https://enzo.readthedocs.org/en/latest/tutorials/building_enzo.html>`_.
To compile Grackle, complete the following procedure:

1. Initialize the build system.

.. highlight:: none

::

    ~ $ cd grackle
    ~/grackle $ ./configure

2. Proceed to the source directory.

.. highlight:: none

::

    ~/grackle $ cd src/clib

3. Configure the build system.

.. note:: 
   As of version 2.1, Grackle uses ``libtool`` for building and installation.  
   As such, both shared and static libraries will be built automatically and 
   it is not necessary to add the -fPIC compiler flag.

Compile settings for different systems are stored in files starting with 
"Make.mach" in the source directory.  Grackle comes with three sample make 
macros: ``Make.mach.darwin`` for Mac OSX, ``Make.mach.linux-gnu`` for 
Linux systems, and an unformatted ``Make.mach.unknown``.  If you have a make 
file prepared for an Enzo install, it cannot be used straight away, but is a 
very good place to start.

Once you have chosen the make file to be used, a few variables should be set:

    * ``MACH_LIBTOOL`` - path to ``libtool`` executable.  Note, on a Mac, 
      this should point to ``glibtool``, which can be installed with macports 
      or homebrew.

    * ``LOCAL_HDF5_INSTALL`` - path to your hdf5 installation.  

    * ``LOCAL_FC_INSTALL`` - path to Fortran compilers (not including the bin 
      subdirectory).

    * ``MACH_INSTALL_PREFIX`` - path where grackle header and library files 
      will be installed.

    * ``MACH_INSTALL_LIB_DIR`` - path where libgrackle will be installed (only 
      set if different from MACH_INSTALL_PREFIX/lib).

    * ``MACH_INSTALL_INCLUDE_DIR`` - path where grackle header files will be 
      installed (only set if different from MACH_INSTALL_PREFIX/include).

Once the proper variables are set, they are loaded into the build system by 
doing the following:

.. highlight:: none

::

    ~/grackle/src/clib $ make machine-<system>

Where system refers to the make file you have chosen.  For example, if you 
chose ``Make.mach.darwin``, type:

.. highlight:: none

::

    ~/grackle/src/clib $ make machine-darwin

Custom make files can be saved and loaded from a **.grackle** directory in the 
home directory.

.. _compiler-settings:

Compiler Settings
+++++++++++++++++

There are three compile options available for setting the precision of 
baryon fields, compiler optimization, and enabling OpenMP.  To see them,
type:

.. highlight:: none

::

    ~/grackle/src/clib $ make show-config

   MACHINE: Darwin (OSX)
   MACHINE-NAME: darwin

   CONFIG_PRECISION  [precision-{32,64}]                     : 64
   CONFIG_OPT  [opt-{warn,debug,high,aggressive}]            : high
   CONFIG_OMP  [omp-{on,off}]                                : off

For example, to change the optimization to high, type:

.. highlight:: none

::

    ~/grackle/src/clib $ make opt-high

.. warning::
   Compiling Grackle in single precision (with ``make precision-32``) is **not**
   recommended. Because of the high dynamic range involved in calculating many
   chemistry and cooling rates, running Grackle in single precision can produce
   unreliable results. This is especially true when running with
   :c:data:`primordial_chemistry` >= 1.

Custom settings can be saved for later use by typing:

.. highlight:: none

::

    ~/grackle/src/clib $ make save-config-<keyword>

They will be saved in the **.grackle** directory in your home directory.  To 
reload them, type:

.. highlight:: none

::

    ~/grackle/src/clib $ make load-config-<keyword>

For a list of all available make commands, type:

.. highlight:: none

::

    ~/grackle/src/clib $ make help

    ========================================================================
       Grackle Makefile Help
    ========================================================================
    
       make                Compile and generate librackle
       make install        Copy the library somewhere
       make help           Display this help information
       make clean          Remove object files, executable, etc.
       make dep            Create make dependencies in DEPEND file
    
       make show-version   Display revision control system branch and revision
       make show-diff      Display local file modifications
    
       make help-config    Display detailed help on configuration make targets
       make show-config    Display the configuration settings
       make show-flags     Display specific compilation flags
       make default        Reset the configuration to the default values

4. Compile and Install

To build the code, type:

::

    ~/grackle/src/clib $ make 
    Updating DEPEND
    Compiling calc_rates.F
    Compiling cool1d_multi.F
    ....
    
    Linking
    Success!

Then, to install:

::

    ~/grackle/src/clib $ make install

5. Test your Installation

Once installed, you can test your installation with the provided example to
assure it is functioning correctly.  If something goes wrong in this process,
check the ``out.compile`` file to see what went wrong during compilation,
or use ``ldd`` (``otool -L`` on Mac) on your executable to determine what went 
wrong during linking.

::

    ~/grackle/src/clib $ cd ../example
    ~/grackle/src/example $ make clean 
    ~/grackle/src/example $ make 

    Compiling cxx_example.C
    Linking
    Success!
  
    ~/grackle/src/example $ ./cxx_example

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

In order to verify that Grackle is fully functional, try :ref:`running the
test suite <testing>`.

.. _cmake_build:

Building with CMake
-------------------

To use this system, version 3.16 or newer of ``cmake`` is required.

.. warning::

   This build-system may not work properly if you have previously tried to build an earlier version of Grackle with the classic build system.

   * While the cmake build system performs an "out-of-source" build, the traditional build system performs an "in-source" build.

   * While the "classic build system" has been modified to better coexist with the cmake build-system, earlier versions could cause issues.
     If the auto-generated files (both headers and source files), produced by the in-source build (from an earlier Grackle-version), are not properly removed, this can cause issues for cmake builds.

Procedure
+++++++++

1. Proceed to the grackle directory

   .. code-block:: shell-session

      ~$ cd grackle


2. Initialize and configure the build-system.
   In these example snippets, we show the minimum required configuration options (this should work on most machines) and provide more details later about :ref:`additional configuration options <available_cmake_options>` and :ref:`how to specify configuration options <how_to_configure>` down below.
   For now, we make 2 basic decisions:

   #. Decide on the directory, ``<build-dir>``, where you want to build Grackle. [#f1]_
      This is referred to as the build-directory and is generally placed at the root level of the grackle repository.
      A common choice is ``build`` (but this is fairly arbitrary).

   #. Decide on the installation directory prefix, ``<install-prefix>``, where Grackle will be installed.
      This is be specified via the ``CMAKE_INSTALL_PREFIX`` cmake configuration variable.
      On UNIX-like systems, it defaults to ``/usr/local/``.

   To configure a build where Grackle is compiled as a static library, use

   .. code-block:: shell-session

      ~/grackle $ cmake -DCMAKE_INSTALL_PREFIX=<install-prefix> -B <build-dir>

   To configure a build where Grackle is compiled as a shared library, use

   .. code-block:: shell-session

      ~/grackle $ cmake -DCMAKE_INSTALL_PREFIX=<install-prefix> -DBUILD_SHARED_LIBS=ON -B <build-dir>

   .. note::

       If you are building Grackle to be used with a downstream simulation-code, that doesn't mention any preferences about how Grackle is built, you will probably have more luck compiling Grackle as a shared library.

3. Compile and install grackle.

   .. code-block:: shell-session

      ~/grackle $ cmake --build <build-dir>
      ~/grackle $ cmake --install <build-dir>

   .. note::

      The above commands show the most generic commands that can be executed.
      Other tutorials that you see online may show slight variations in these commands (where you manually make the build directory) and then manually execute the build-system from within the build-directory...

   .. note::

      Just like with the classic build-system, Grackle currently needs to be installed to be used.
      If you install it in a non-standard location, then you also need to ensure that you properly set the LD_LIBRARY_PATH (or DYLD_LIBRARY_PATH on macOS) to make use of it.

      The current structure (and contents) of the build-directory can and will change (especially once `GH-#204 <https://github.com/grackle-project/grackle/pull/204>`__ and `GH-#208 <https://github.com/grackle-project/grackle/pull/208>`__ are merged).
      But, we plan to add support for linking against a grackle installation without fully installing it.


4. Test your Build.

   Once you have compiled Grackle, you can run one of the provided example to test if it functions correctly.
   These examples are automatically compiled with Grackle.

   .. code-block:: shell-session

      ~/grackle $ cd <build-dir>/examples
      ~/grackle/<build-dir>/examples $ ./cxx_example

   .. warning::

      The examples make certain assumptions about the location of the input files.
      The examples are only guarunteed to work if both:

         1. you execute the example-binary from the same-directory where the example-binary is found

         2. ``<build-dir>`` is a top-level directory in the grackle repository (e.g. something like ``my-build`` is fine, but choices like ``../my-grackle-build`` and ``my_builds/my-first-build`` are problematic).

   .. note::

      For reference, the Classic build-system always links Grackle against the shared-library version of Grackle and requires that Grackle is fully installed in a location known by the system (either a standard system location OR a location specified by ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH``).
      In contrast, cmake automatically takes special-steps to try to ensure that each example-binary will link to the copy of the Grackle library (whether it is shared or static) that is in the ``<build-dir>``; in fact, Grackle doesn't even need to be installed to run the Grackle library.

      With that said, if you compile Grackle as a shared library in a cmake build, an example-binary **might** try to use a copy of a shared grackle library found in a directory specified by ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH`` if one exists.
      The exact behavior may be platform dependent and also depends on whether CMake instructs the linker to use RPATH or RUNPATH (this is not spacified by the cmake docs).

.. _how_to_configure:

How to Specify Configuration Options
++++++++++++++++++++++++++++++++++++

All configuration options can be specified when invoking cmake during configuration of the build.
Specifically you can specify can specify the values by inserting an argument of the form ``-D<variable>=<value>`` to the list of arguments passed to ``cmake``.
This is illustrated in the prior subsection where we pass ``-DCMAKE_INSTALL_PREFIX=/my/install/path...`` and ``-DBUILD_SHARED_LIBS=OFF``.

Alternatively, you can replace the call to ``cmake`` during configuration with a call to ``ccmake`` to provide a TUI (text-based user interface) where you can manually configure options.
For example, a call to ``ccmake -B<build-dir>`` will bring up a TUI to configure a build in the specified directory.
CMake also provides a GUI (graphical user interface) for this purpose (it may not be available based on how exactly you installed CMake).
The CMake documentation provide more details about the GUI and how to more generally use cmake `here <https://cmake.org/cmake/help/latest/guide/user-interaction/index.html#guide:User%20Interaction%20Guide>`__.

A summary of all Grackle-specific configuration options and a subset of useful generic CMake configurations is provided in the :ref:`next subsection <available_cmake_options>`.

The idiomatic way to control optimization/debugger flags is store a build-type in the standard ``CMAKE_BUILD_TYPE`` variable.
Choices include:

* ``-DCMAKE_BUILD_TYPE=Release`` (typically ``-O3``)

* ``-DCMAKE_BUILD_TYPE=RelWithDebInfo`` (typically ``-O2 -g``)

* ``-DCMAKE_BUILD_TYPE=Debug`` (typically ``-O0 -g``)

The first choice is generally fastest, while the second is a sensible choice during development (the compiler performs most optimizations and includes debugging information in the library).


*[ NEED TO ADDRESS: machine files and* ``CMAKE_<LANG>_FLAGS`` *]*

.. _available_cmake_options:

Available Configuration Options
+++++++++++++++++++++++++++++++

The compilation (and installation) of Grackle can be configured using various options.
These options are described in the following 2 tables.

Many of these options are binary choices that accept a boolean value. [#f2]_

This first table describes the Grackle-specific options to configure the build.

.. list-table:: Grackle-Specific Options
   :widths: 12 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``GRACKLE_USE_DOUBLE``
     - Turn off to build Grackle with single precision.
     - ``"ON"``
   * - ``GRACKLE_USE_OPENMP``
     - Turn on to build Grackle with OpenMP
     - ``"OFF"``

This second table highlights a subset of standardized CMake options that may also be useful.

.. list-table:: Standard CMake Options
   :widths: 12 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default

   * - ``BUILD_SHARED_LIBS``
     - When turned ``"ON"``, Grackle is built as a shared library. When turned ``"OFF"`` (or if its undefined), Grackle is built as a static library.
     - ``<undefined>``

   * - ``CMAKE_BUILD_TYPE``
     - Specifies the desired build configuration (for single-configuration generators [#f3]_).
       Grackle currently supports the standard choices ``Debug``, ``Release``, ``RelWithDebInfo`` and ``MinSizeRel``.
     - ``<undefined>``

   * - ``CMAKE_INSTALL_PREFIX``
     - Specifies the path-prefix where Grackle will be installed when you invoke ``make install`` from within the build-directory (or using a non-Makefile generator, you use the generator-specific command to build the ``install``-target).
       Note, that if you use ``cmake --install path/to/builddir`` to invoke installation, you can use ``--prefix`` to specify a different prefix
     - ``/usr/local``

   * - ``HDF5_ROOT``
     - When cmake has trouble finding your hdf5 installation, you can set this variable equal to the path to the HDF5 installation to serve as a hint for cmake
     - ``<undefined>``

   * - ``HDF5_PREFER_PARALLEL``
     - Set to ``true`` to express a preference for linking against parallel hdf5 (by default, the serial version will be preferentially choosen)
     - ``<undefined>``

There are also additional standard options for BOTH configuring other aspects of the build and for finding the correct/preferred HDF5 library and configuring the correct openmp library.



.. _cmake_shared_and_static:

Installing both Shared and Static Libraries
+++++++++++++++++++++++++++++++++++++++++++

It's idiomatic for a given ``cmake``-build to build either a shared library OR a static library (not both). This is controlled by the standard ``BUILD_SHARED_LIBS`` flag (you usually don't need both).

With that said, if you really want to install both of them, you could trigger 2 separate builds that install to the same destination. [#f4]_
The following code snippet illustrates how you might do this (for concreteness, the snippet uses build-directories called ``build-static`` and ``build-shared`` and installs into a directory called ``$HOME/local`` -- but these are all arbitrary choices).
   
.. code-block:: shell-session

   ~ grackle $ cmake -DCMAKE_INSTALL_PREFIX=$HOME/local -B build-static
   ~ grackle $ cmake --build build-static
   ~ grackle $ cmake --install build-static
   ~ grackle $ cmake -DCMAKE_INSTALL_PREFIX=$HOME/local -DBUILD_SHARED_LIBS=ON -B build-shared
   ~ grackle $ cmake --build build-shared
   ~ grackle $ cmake --install build-shared

Purpose
+++++++

The purpose of this build-system is to facillitate more seamless integration with downstream applications built with CMake.
In particular, our primary focus is to allow developers to directly embed Grackle into their application.
This is commonly achieved with git submodules.

There are a few benefits to this approach:

- This integration makes for a very streamlined installation experience for end-users.

  - To install the downstream application, the end-user just has to:

      1. clone the downstream application and initialize all git submodules

      2. configure and build the downstream application.

    The downstream application is able to automatically include the compilation of Grackle as part of its build process.

  - The above process includes far fewer steps than the more traditional installation process.
    In the more traditional procedure, the user must (i) clone Grackle and (ii) configure and build Grackle before they execute the steps in the above bullet.
    They must also worry about configuring the installation of the downstream application to properly find Grackle.

- This approach also simplifies scripts used for automated testing of the downstream.
  Obviously, the streamlined installation process will simplify the scripts (especially if you want to try compiling with single vs. double precision).
  There is also a more subtle benefit: the grackle datafiles have a predictable location.

Two considerations that should be weighed before considering this approach:

1. The downstream application's build system needs to include some extra logic to properly configure the Grackle build.
   In reality, many/all of Grackle's dependencies are probably already dependencies of the downstream application (e.g. hdf5).
   Additionally, this logic of say choosing Grackle's floating-point precision may be able to replace existing compatability checks.

2. Developers of downstream application need to keep updating the `.gitmodules` file as newer grackle versions are released.
   If the developers are already following best practices, this probably isn't much extra work.
   Ideally, they should already be informing their users about grackle version compatability.
   The `.gitmodules` file can be considered a centralized location where this compatability can be checked.

.. COMMENT-BLOCK

   **As an aside:** this embedding approach will implicitly encourage downstream developers to avoid the common pitfall in CI scripts of simply downloading the most recent release of a dependency.

Finally, it's worth mentioning that a downstream project can be configured to use either this-embedded approach **OR** link against a separately compiled version of Grackle (using the more traditional build-system).
This is probably advisable, given the experimental nature of this buildsystem (e.g. so if a user runs into problems with the CMake-build of Grackle on some more uncommon system, they can always fall back to the more traditional approach).
We will discuss how to do it down below.

.. note::

   This remainder of this section assumes that the reader already has familiarity with ``cmake``.



.. rubric:: Footnotes

.. [#f1] For the uninitiated, Grackle performs "out of source builds," in which the build-artifacts, like generated headers, object files, linked libraries, are placed inside a build directory (rather than putting them inside the source-directory next to the source files).
         There are a couple of advantages to this approach such as (i) you can maintain multiple builds at the same time (e.g. if you are switching between development branches) or (ii) it's really easy to clean up from a build (you just delete the build-directory).



.. [#f2] CMake boolean variables map a variety of values to ``true`` (e.g. ``1``, ``ON``, ``TRUE``, ``YES``, ``Y``) and a variety of values to ``false`` (e.g. ``0``, ``OFF``, ``FALSE``, ``NO``, ``N``).

.. [#f3] If you are simply following the above compilation instructions, you definitely don't need to worry about the distinction between a single-configuration generator (e.g. Makefiles and standard Ninja) and multi-configuration generators.

.. [#f4] Aside: performing these 2 separate CMake builds compiles the source files the same number of times as the Classic build system.
         Behind the scenes, the classic build system always compile each source file twice (once with position independent code and once without).



