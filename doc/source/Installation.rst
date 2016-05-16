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

Grackle is available in a mercurial repository 
`here <https://bitbucket.org/grackle/grackle>`__.  The mercurial site 
is `here <http://mercurial.selenic.com/>`__ and an excellent tutorial can be 
found `here <http://hginit.com/>`__.  With mercurial 
installed, grackle can be obtained with the following command:

.. highlight:: none

::

    ~ $ hg clone https://bitbucket.org/grackle/grackle

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
