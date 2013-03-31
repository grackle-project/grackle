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
`here <https://bitbucket.org/brittonsmith/grackle>`_.  The mercurial site 
is `here <http://mercurial.selenic.com/>`_ and an excellent tutorial can be 
found `here <http://hginit.com/>`_.  With mercurial 
installed, grackle can be obtained with the following command:

.. highlight:: none

::

    ~ $ hg clone https://bitbucket.org/brittonsmith/grackle

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

Compile settings for different systems are stored in files starting with 
"Make.mach" in the source directory.  Grackle comes with three sample make 
macros: ``Make.mach.darwin`` for Mac OSX, ``Make.mach.linux-gnu`` for 
Linux systems, and an unformatted ``Make.mach.unknown``.  If you have a make 
file prepared for an Enzo install, you may use it to compile grackle.
Once you have chosen the make file to be used, a few variables should be set:

    * ``LOCAL_HDF5_INSTALL`` - path to your hdf5 installation.

    * ``LOCAL_FC_INSTALL`` - path to fortran compilers (not including the bin subdirectory).

    * ``MACH_INSTALL_PREFIX`` - path where grackle header and library files will be installed.

    * ``MACH_INSTALL_LIB_DIR`` - path where libgrackle will be installed (only set if different from MACH_INSTALL_PREFIX/lib).

    * ``MACH_INSTALL_INCLUDE_DIR`` - path where grackle header files will be installed (only set if different from MACH_INSTALL_PREFIX/include).

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

Compiler Settings
+++++++++++++++++

There are three compile options available for setting the precision of 
floating point and integer variables and for optimization.  To see them, 
type:

.. highlight:: none

::

    ~/grackle/src/clib $ make show-config

   MACHINE: Darwin (OSX)
   MACHINE-NAME: darwin

   CONFIG_PRECISION  [precision-{32,64}]                     : 64
   CONFIG_INTEGERS  [integers-{32,64}]                       : 64
   CONFIG_OPT  [opt-{warn,debug,cudadebug,high,aggressive}]  : debug

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

Now it's time to integrate grackle into your simulation code: :ref:`integration`
