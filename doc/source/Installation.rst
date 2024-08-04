.. _obtaining_and_building_enzo:

Installation
============

There are 3 steps to setting up Grackle on your system

   1. :ref:`Install Grackle's Dependencies <install_grackle_dependencies>`

   2. :ref:`Download Grackle <download_grackle>`

   3. Build and install Grackle using the :ref:`classic build system <classic_build>` or the :ref:`CMake build system <cmake_build>`.


.. note::

   Given a smooth roll-out of the :ref:`CMake build system <cmake_build>`, it is our intention to deprecate and remove the :ref:`classic build system <classic_build>`.
   If you encounter any problems with the CMake system or anticipate any issues with this plan, :doc:`please let us know <Help>`.

We include a :ref:`note on compiler toolchain compatability <compiler_toolchain_compatability>` at the end of this page.


.. _install_grackle_dependencies:

Dependencies
------------

In addition to C/C++ and Fortran compilers, the following dependency must 
also be installed:

   * `HDF5 <http://www.hdfgroup.org/HDF5/>`_, the hierarchical data format.
     HDF5 also may require the szip and zlib libraries, which can be
     found at the HDF5 website.

     * For the :ref:`classic build system <classic_build>`, compiling with HDF5 1.8 or greater requires that the ``H5_USE_16_API`` compatability directive is manually specified.
       This can be done by adding ``-DH5_USE_16_API`` to the list of compiler flags given in machine-specific make files.

     * The :ref:`CMake build system <cmake_build>`, automatically handles these details for you.

Although many systems already have them installed, both build systems have additional dependencies:

   * the :ref:`classic build system <classic_build>`, employs the ``makedepend`` and the `libtool <https://www.gnu.org/software/libtool/>`_ utilities.
     It's often easiest to download these dependencies through your system's package manager.

   * the :ref:`CMake build system <cmake_build>` requires cmake to be installed.
     It's easiest download a binary distribution from the `CMake website <https://cmake.org/download/>`_ or use your system's package manager.
     We require version 3.16 or newer.

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

Grackle provides a Modern CMake build-system.
While CMake has some baggage (primarily due to the maintenace of backwards compatability), it is arguably the most-portable mainstream build-system that is easiest to integrate with simulation codes.

An overview of our design philosophy is provided :ref:`here <cmake_buildsystem_design_rationale>`.
This build-system makes integration of Grackle into simulation codes that are themselves built with CMake extremely easy.
Steps have also been taken simplify integration of Grackle into simulation codes built with any other build-systems (they just need to call the standardized ``pkg-config`` command-line tool).
More details about integration are provided :doc:`on this page <Consuming>`.
This current section focuses on installation.

For the uninitiated, the CMake build-system performs an out-of-source build.
An out-of-source build places all build artifacts (auto-generated source/header files, object files, etc.) into a "build-directory."
The build-directory is at a user-specified location that is organized into a hierarchy that resembles the source directory hierarchy.
Cleaning up from a CMake-build is as simple as deleting this build-directory.
In contrast, the "classic build system" performs an in-source build (because that type of build distributes build artifacts throughout the source directory hierarchy, clean up requires more complex logic encapsulated by the ``make clean`` command).

.. warning::

   While the "classic build system" has been modified to better coexist with the CMake build-system, issues can potentially arise if build-artifacts produced in a "classic" build of an earlier Grackle-revision are not properly removed.
   Specifically, the issues relate to the presence of auto-generated header-files.
   We have built checks into the CMake build-system to prevent these issues in most cases, but they may not help in certain pathological scenarios.

Procedure
+++++++++

The build/installation procedure follows the standard steps of any CMake build.
The remainder of this subsection is primarily intended for readers who are relatively inexperienced with using CMake.

1. Proceed to the grackle directory

   .. code-block:: shell-session

      ~$ cd grackle


2. Initialize and configure the build-system.
   In these example snippets, we show the minimum required configuration options (this should work on most machines) and provide more details later about :ref:`additional configuration options <available_cmake_options>` and :ref:`how to specify configuration options <how_to_configure>` down below.
   During this step you might also specifiy :ref:`machine-specific host files <cmake_host-files>` (but that usually isn't absolutely necessary).

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


   It is idiomatic for a given CMake build to just compile Grackle as either a static or shared library, not both (you usually just need one).
   But if you must have both, see :ref:`this section <cmake_shared_and_static>`.

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

      The current structure (and contents) of the build-directory can and will change (especially once `GH-#208 <https://github.com/grackle-project/grackle/pull/208>`__ is merged).
      But, we plan to add support for linking against a grackle installation without fully installing it.


4. Test your Build.

   Once you have compiled Grackle, you can run one of the provided example to test if it functions correctly.
   These examples are automatically compiled with Grackle.

   .. code-block:: shell-session

      ~/grackle $ cd <build-dir>/examples
      ~/grackle/<build-dir>/examples $ ./cxx_example

   .. warning::

      The examples make certain assumptions about the location of the input files.
      The examples are only guaranteed to work if both:

         1. you execute the example-binary from the same-directory where the example-binary is found

         2. ``<build-dir>`` is a top-level directory in the grackle repository (e.g. something like ``my-build`` is fine, but choices like ``../my-grackle-build`` and ``my_builds/my-first-build`` are problematic).

   .. note::

      For reference, the Classic build-system always links Grackle against the shared-library version of Grackle and requires that Grackle is fully installed in a location known by the system (either a standard system location OR a location specified by ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH``).
      In contrast, cmake automatically takes special-steps to try to ensure that each example-binary will link to the copy of the Grackle library (whether it is shared or static) that is in the ``<build-dir>``; in fact, Grackle doesn't even need to be installed to run the Grackle library.

      With that said, if you compile Grackle as a shared library in a cmake build, an example-binary **might** try to use a copy of a shared grackle library found in a directory specified by ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH`` if one exists.
      The exact behavior may be platform dependent and also depends on whether CMake instructs the linker to use RPATH or RUNPATH (this is not specified by the cmake docs).

In order to verify that Grackle is fully functional, you can try :ref:`running the test suite <testing>`.

.. _how_to_configure:

How to Specify Configuration Options
++++++++++++++++++++++++++++++++++++

All configuration options can be specified when invoking cmake during configuration of the build.
Specifically you can specify the values by inserting an argument of the form ``-D<variable>=<value>`` to the list of arguments passed to ``cmake``.
This is illustrated in the prior subsection where we pass ``-DCMAKE_INSTALL_PREFIX=/my/install/path...`` and ``-DBUILD_SHARED_LIBS=OFF``.

Alternatively, you can replace the call to ``cmake`` during configuration with a call to ``ccmake`` to provide a TUI (text-based user interface) where you can manually configure options.
For example, a call to ``ccmake -B<build-dir>`` will bring up a TUI to configure a build in the specified directory.
CMake also provides a GUI (graphical user interface) for this purpose (it may not be available based on how exactly you installed CMake).
The CMake documentation provide more details about the GUI and how to more generally use cmake `here <https://cmake.org/cmake/help/latest/guide/user-interaction/index.html#guide:User%20Interaction%20Guide>`__.

A summary of all Grackle-specific configuration options and a subset of useful generic CMake configurations is provided in the :ref:`next subsection <available_cmake_options>`.

The idiomatic way to control optimization/debugger flags is to specify a build-type via the standard ``CMAKE_BUILD_TYPE`` variable.
Choices include:

* ``-DCMAKE_BUILD_TYPE=Release`` (typically ``-O3``)

* ``-DCMAKE_BUILD_TYPE=RelWithDebInfo`` (typically ``-O2 -g``)

* ``-DCMAKE_BUILD_TYPE=Debug`` (typically ``-O0 -g``)

The first choice is generally fastest, while the second is a sensible choice during development (the compiler performs most optimizations and includes debugging information in the library).

Machine-specific compilation options can also be specified with host-files.
These host-files should generally not be necessary, but they may specify architecture-specific optimization flags.
This should be specified during the configuration stage with the ``-C`` flag followed by the path to the host-file.
For example, one might invoke:

   .. code-block:: shell-session

      ~/grackle $ cmake -C config/host-config/tacc-frontera-intel.cmake \
      > -D CMAKE_INSTALL_PREFIX=<install-prefix> \
      > -D BUILD_SHARED_LIBS=ON \
      > -B <build-dir>

The order of ``-D`` and ``-C`` flags matters.
If they are both used to specify values for a given variable, the last one to appear "wins."
More information about writing host-files are provided :ref:`below <cmake_host-files>`.


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

   * - ``CMAKE_<LANG>_COMPILER``
     - Set of variables (where ``<LANG>`` is replaced by ``C``, ``Fortran`` or ``CXX``) to overide the compiler choice.
       This is commonly set by host-files.
     - ``<undefined>``

There are also additional standard options for BOTH configuring other aspects of the build and for finding the correct/preferred HDF5 library and configuring the correct openmp library.

Addtionally, CMake will also respect the values of certain environment variables.
For example, if you don't manually specify the choice of compilers with the ``CMAKE_<LANG>_COMPILER`` flag, then CMake will use the values in the ``CC``, ``FC``, and ``CXX`` environment variables.

We strongly encourage users and developers to make use of the options described in this section.
They exist to provide a curated/consistent experience in a variety of scenarios.
:doc:`Please let us know <Help>` if you think we are missing a useful Grackle-specific option.
You can also add the new option yourself (it may be useful to review :ref:`the design philosophy for the CMake build-system <cmake_buildsystem_design_rationale>`).

With that said, we also recognize that the need may arise where a user/developer may want to specify arbitrary flags.
You can use the standardized ``CMAKE_<LANG>_FLAGS`` variables for that purpose (where ``<LANG>`` is ``C``, ``CXX``, ``Fortran``).
For example, passing ``-DCMAKE_C_FLAGS="-Wall -Wpedantic -funroll-loops"`` will pass these flags to every invocation of the C compiler (for compiling Grackle itself as well as any examples or tests).
Technically, these flags are passed to every invocation of the C compiler-frontend (even during linking), but that usually isn't a problem.



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

.. _cmake_host-files:

More About Host-Files
+++++++++++++++++++++

As noted above, we provide support for setting default value for particular machines by providing support for *host-files*\ .
These files are provided mostly for convenience (and to provide parity with machine files provided by the classic build-system).
They are most useful on HPC systems that provide multiple compiler toolchains.
These are the *\*.cmake* files in the **config/host-config** directory.

Importantly, the usage of *host-files* is optional (and usually not required).
They usually aren't needed on local systems (if you find that Grackle won't compile without a host-file, please let us know -- that may indicative of a bug).
They should generally **NOT** be used when Grackle is embedded within another CMake project.

While there are a couple of ways to implement this concept, our current strategy draws some inspiration from `here <https://llnl-blt.readthedocs.io/en/develop/tutorial/host_configs.html>`__.
Essentially, our strategy leverages cmake functionality to pre-load a script to populate some cache variables.

Usually, will specify the desired compilers.
If a HPC machine properly manages the ``CC``, ``FC``, and ``CXX`` environment variables this isn't strictly necessary.
If the machine places HDF5 in an unusual location, you might also hardcode hints into the config-file.

The most important role is to specify cluster-specific optimization flags via the special Grackle-specific ``GRACKLE_OPTIMIZATION_FLIST_INIT`` variable.
These flags will **ONLY** be used when compiling Grackle with the ``Release`` or ``RelWithDebInfo`` build-types.
Here are 2 illustrative examples:

 * First we show that in order to pass multiple flags, the flags need to be specified by a semicolon delimited list.
   If you stored ``"-xCORE-AVX512;-funroll-loops"`` within ``GRACKLE_OPTIMIZATION_FLIST_INIT``, then all source files will be compiled with these options (they won't be passed to the linker).

 * Next we show that to properly pass "option groups" you may need to make use of CMake's shell-like quoting with the ``SHELL:`` prefix (this relates to option de-duplication performed by CMake).
   Thus, storing ``"SHELL:-option1 A;-Wall;SHELL:-option2 B"`` within ``GRACKLE_OPTIMIZATION_FLIST_INIT`` would cause all compiler invocations for source files used in Grackle to be passed ``-option1 A -Wall -option2 B``.

While embedded builds currently respect ``GRACKLE_OPTIMIZATION_FLIST_INIT``, that is something we may stop supporting.

.. COMMENT-BLOCK

   The tone of this section should make it clear that host-files usually aren't necessary in most scenarios.
   There's a chance that may change if we start supporting CUDA or HIP, these may become more important.
   Until then, I'm a little hesitant to really encourage them since it may unnecessarily complicate things.

.. note::

   If you want to pass language-specific optimization options, let us know.
   That is something we can easily support.
   Until then, this could be addressed by enclosing a given option (or option-group) within a `language-specific generator expressions <https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html#genex:COMPILE_LANGUAGE>`__.

.. note::

   In terms of modern, idiomatic CMake features, a host-file could be replaced by a combination of a `toolchain-file <https://cmake.org/cmake/help/latest/manual/cmake-toolchains.7.html>`__ and a `preset-file <https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html>`__.

   * toolchain files usually define compiler-toolchain related-options and are commonly used for cross-compiling. 
     As a basic rule of thumb: you should be able to recycle toolchain-files between unrelated projects (i.e. they don't include project-specific variables)

   * a preset file (``CMakePresets.json`` or ``CMakeUserPresets.json``) is intended to be used to specify common project-specific compilation options.
     These can be read by IDEs.

   * after we update the minimum required CMake version for compiling Grackle to at least 3.19, we may transition to using these features.


.. _compiler_toolchain_compatability:

Compiler Toolchain Compatability
--------------------------------

As a general rule of thumb, the easiest, most reliable thing to do is  to ensure that Grackle is built with the same compiler toolchain (or a compatible one) as the

   * the downstream application itself (whether it's a simulation code or pygrackle)
   * any other dependencies of the application (whether it's other software libraries or other python extension-modules loaded at the same time).

This is only something you need to consider on platforms with multiple compiler toolchains present. 

In practice, toolchain-compatibility generally **ISN'T** much of a concern for Grackle, when compiled without OpenMP.
In this scenario, you need to use Fortran compilers with consistent runtime libraries (e.g. you might encounter issues if you use ``gfortran`` to compile Grackle and ``ifort`` to compile a downstream simulation code).
If the downstream application doesn't use any Fortran, then there generally aren't any concerns at all.

Things are slightly more complex when compiling Grackle with OpenMP.
You need to make sure that your C compiler and Fortran compiler use a compatible OpenMP runtime.
Usually, your best bet is to try to use C and Fortran compilers from the same vendor (e.g. using ``gcc`` with ``gfortran`` will work or using ``icc`` with ``ifort`` will work).
You might be able to mix compilers from different vendors by passing special compiler and linker options, but this usually isn't well documented.
If your downstream application is also compiled with OpenMP, you also need to ensure that the downstream application is compiled with a compatible runtime.

You don't generally need to worry about OpenMP-compatability between Grackle and the rest of the software stack if Grackle is compiled without OpenMP or if it is the only part of the software stack that is compiled with OpenMP.

**As Grackle continues to evolve, compiler toolchain compatability will become more of an issue.**
For example, adding GPU-support with the likes of CUDA or HIP would involve linking to a C++ runtime library.

.. note::

   Mixing compiler toolchains may be more difficult for certain vendors.
   For example, some vendors may more aggressively link their OpenMP runtime library or C++ runtime-library libraries to the resulting binaries, which could easily cause problems.
   But generally, GNU-compilers and clang are pretty good about this.


.. rubric:: Footnotes

.. [#f1] For the uninitiated, Grackle performs "out of source builds," in which the build-artifacts, like generated headers, object files, linked libraries, are placed inside a build directory (rather than putting them inside the source-directory next to the source files).
         There are a couple of advantages to this approach such as (i) you can maintain multiple builds at the same time (e.g. if you are switching between development branches) or (ii) it's really easy to clean up from a build (you just delete the build-directory).



.. [#f2] CMake boolean variables map a variety of values to ``true`` (e.g. ``1``, ``ON``, ``TRUE``, ``YES``, ``Y``) and a variety of values to ``false`` (e.g. ``0``, ``OFF``, ``FALSE``, ``NO``, ``N``).

.. [#f3] If you are simply following the above compilation instructions, you definitely don't need to worry about the distinction between a single-configuration generator (e.g. Makefiles and standard Ninja) and multi-configuration generators.

.. [#f4] Aside: performing these 2 separate CMake builds compiles the source files the same number of times as the Classic build system.
         Behind the scenes, the classic build system always compile each source file twice (once with position independent code and once without).



