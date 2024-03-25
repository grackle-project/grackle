.. _obtaining_and_building_enzo:

Installation
============

The compilation process for grackle is very similar to that for 
`Enzo <http://enzo-project.org>`_.  For more details on the Enzo build 
system, see the `Enzo build documentation 
<https://enzo.readthedocs.org/en/latest/tutorials/building_enzo.html>`_.  

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

Building
--------

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

Experimental CMake Build System
-------------------------------

Recent versions of Grackle include a secondary, experimental build-system depending on cmake (it's intended to supplement the primary/traditional build system in special cases).
To use this system, version 3.16 or newer of ``cmake`` is required (**NOTE:** ``cmake`` is **NOT** required if you are using the primary build-system).

.. warning::

   This build-system may not work properly if you have previously tried to build grackle with the traditional build system.
   While the cmake build system performs an "out-of-source" build, the traditional build system performs an "in-source" build.
   Specifically, some of the auto-generated files (both headers and source files), produced by the in-source build, can cause issues for cmake builds.

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

Configuring Grackle's Build
+++++++++++++++++++++++++++

The following table lists Grackle-specific cmake options that can be used to configure the build

.. list-table:: General Configuration
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``GRACKLE_USE_DOUBLE``
     - Turn off to build Grackle with single precision.
     - ON
   * - ``GRACKLE_USE_OPENMP``
     - Turn on to build Grackle with OpenMP
     - OFF

You can use standard options to provide the build system with hints for finding the correct HDF5 library and configuring the correct openmp library.

There are 2 noteworthy differences from the traditional build system:

1. It's idiomatic for a given ``cmake`` to build either a shared library OR a static library (not both). This is controlled by the standard ``BUILD_SHARED_LIBS`` flag.

2. (On at least some platforms) When ``cmake`` constructs a shared libraries with ``OPENMP`` support, the resulting library is "more fully" linked against the OPENMP runtime library.
   Downstream applications don't need to know anything about whether such a Grackle library uses OpenMP during compilation (this contrasts with the more traditional approach, where you would explicitly need to link against openmp).
   This has an interesting consequence that you could compile pygrackle with openmp support.

Instructions for Integration
++++++++++++++++++++++++++++

*[ TO BE ADDED ]*

.. COMMENT-BLOCK

   I need to actually test this all out. The crux of this is add_subdirectory and then link against Grackle::Grackle
   I'll definitely circle back and update this.

.. note::

   At this time, we do NOT support ordinary installation with the cmake build-system.
   Setting this up right (especially while promoting compatibility with the embedding approach) is a little challenging.
   The CMake provided machinery for this task is somewhat complex, because the task of supporting installation various kinds of installations across multiple platforms is itself complex.
