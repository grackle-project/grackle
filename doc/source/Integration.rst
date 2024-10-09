Integrating Grackle into your Application
=========================================

**If you are simply installing Grackle because it is a dependency of a simulation-code that you want to install (and you aren't a developer of that code), you can probably disregard this section.**

.. note::

   This page employs some strong language to express the developers' intentions about what are *officially* supported ways for Grackle to be used.
   Be advised that things may break if you use Grackle in an unsupported manner.

   If there is some other way you want to consume Grackle (that isn't listed here), please ping us and let us know!
   We are happy to try to accommodate your requests!
   But please try to let us know sooner rather than later.
   If we haven't anticipated the particular way that you want to consume the library, we may want to make few quick, minor tweaks to make it easier for us to smoothly support your preference as Grackle continues to evolve.

This section assumes that the reader wants to integrate Grackle into their simulation code.
The following table summarizes the *officially* supported approaches that your simulation code could take to support Grackle.
Each of these approaches is compatable with using Grackle's CMake build-system.


+----------------------------------+-----------------------------+------------------+--------------------------+------------------------+
| Name                             | Required downstream         | Supports         | | Works                  | Other Notes            |
|                                  | build-system                | :ref:`classic    | | without                |                        |
|                                  |                             | <classic_build>` | | full install           |                        |
|                                  |                             | builds           |                          |                        |
+==================================+=============================+==================+==========================+========================+
| :ref:`link: manual               | ANY                         | YES              | NO                       | We don't *officially*  |
| <manual_grackle_linking>`        |                             |                  |                          | support this case with |
|                                  |                             |                  |                          | static Grackle libs    |
+----------------------------------+-----------------------------+------------------+--------------------------+------------------------+
| :ref:`link: pkg-config           | ANY\*                       | NO               | NO                       | This is the            |
| <pkgconfig_grackle_linking>`     |                             |                  |                          | recommended way to     |
|                                  |                             |                  |                          | consume Grackle if     |
|                                  |                             |                  |                          | not using CMake        |
+----------------------------------+-----------------------------+------------------+--------------------------+------------------------+
| :ref:`link: CMake                | CMake                       | NO               | EXPERIMENTAL             | We generally recommend |
| <cmake_grackle_linking>`         |                             |                  |                          | the next option        |
|                                  |                             |                  |                          | over this approach     |
+----------------------------------+-----------------------------+------------------+--------------------------+------------------------+
| :ref:`Embedded build             | CMake                       | This is a special case! The end-user doesn't directly build Grackle. |  
| <Embed_Grackle_in_Sim_Build>`    |                             | Instead, Grackle is built as part of the downstream project.         |
|                                  | (`Meson?`_)                 |                                                                      |
|                                  |                             | **This is the most convenient choice for end-users.**                |
+----------------------------------+-----------------------------+------------------+--------------------------+------------------------+

.. _Meson?: https://mesonbuild.com/CMake-module.html#cmake-subprojects


Requirements for Integration Approaches
---------------------------------------

We encourage developers of downstream applications to clearly document how they expect their end-users to install Grackle.

If you want to support end-users that are installing Grackle with the classic build-system, there is currently only 1 *officially* :ref:`supported way <manual_grackle_linking>` to do this.
This approach also supports the case where the end-user installed Grackle as a shared library with the CMake build-system.

All of the other approaches involve compiling Grackle with the CMake build-system and provide *official* support for consuming Grackle as a static library. 

.. note::

   In a vacuum, it is generally preferable to consume Grackle as a static library because the resulting code will be marginally faster, you/your users don't need to worry about managing ``-rpath`` linker-flags or the LD_LIBRARY_PATH variable when Grackle isn't installed in a standard system location, and the Grackle library can be freely moved/deleted after compiling your simulation without breaking your simulation.
   Additionally, the current structure of Grackle's API creates a situation where you should always recompile your simulation whenever Grackle is modified, even if you are consuming Grackle as a shared library (eliminating a strength of shared libraries).

   The only advantage to consuming Grackle as a shared library is that Grackle's implicit runtime dependencies are more easily linked.
   While linking issues should **NEVER** occur for the :ref:`embedding approach <Embed_Grackle_in_Sim_Build>`, there is a small chance they could come up for the remaining approaches (please open an issue if that happens!).
   In the event linking issues do arise, the remaining approaches make it is easy to seamlessly switch from consuming Grackle as a static-library to consuming it a shared library (the end-user can simply delete the existing installation and reinstall it as a shared library).


.. _manual_grackle_linking:

Manually Specifying Linking Flags (when using Grackle as a shared library)
--------------------------------------------------------------------------

This is the historic approach that we have always supported.
If you employ this approach, you should inform your end-users that they should employ, build, and install Grackle via the classic build system or use the CMake build-system to build (and install) Grackle as a shared library.

Here is a minimal sample Makefile for compiling the ``c_example.c`` file.

.. code-block:: makefile

   # To support the case where grackle installations in non-standard locations
   # (without requiring LD_LIBRARY_PATH), you could add 
   # -Wl,-rpath,${GRACKLE_INSTALL}/lib to the start of the LDFLAGS

   # Before Grackle 3.2, CFLAGS needed to include -DCONFIG_BFLOAT_4 or
   # -DCONFIG_BFLOAT_8 based on how grackle was compiled

   GRACKLE_INSTALL = /path/to/grackle/install
   CFLAGS = -I${GRACKLE_INSTALL}/include
   LDFLAGS = -L${GRACKLE_INSTALL}/lib -lgrackle

   # flags unrelated to Grackle
   UNAME := $(shell uname)
   ifneq ($(UNAME), Darwin)
     OTHER_LDFLAGS=-lm
   endif

   c_example:
           $(CC) $(CFLAGS) -c c_example.c -o c_example.o
           $(CC) $(LDFLAGS) $(OTHER_LDFLAGS) c_example.o -o c_example

This scenario has simple linker flags because we are linking to the shared library form of Grackle.
If you used Grackle in its static library form, you would also need to manually link to all of Grackle's dependencies (or at least the ones that don't overlap with your simulations other dependencies.
Some of these dependencies are implicit and depend on the precise choice of compiler and whether openmp is used.

.. warning::

   If a Grackle installation using the Classic Build system made use of an incomplete machine file, the Grackle shared library may not be properly linked to all required runtime dependencies.
   Consequently, the simulation code may need to link against extra dependencies
   This commonly happens with the OpenMP runtime libraries.

   The CMake-builds are much more robust against these kinds of errors.

.. note::

   While you are welcomed to try to link Grackle as a static-library, this is **not** an approach we can **officially** support.
   Complications arise because the set of extra dependencies that must be manually linked is platform-dependent.
   For simulation codes, built with a non-CMake build-system, we recommend :ref:`the pkg-config approach <pkgconfig_grackle_linking>`.

.. _pkgconfig_grackle_linking:

Using pkg-config
----------------

.. note::

   This approach **ONLY** works if the end-user built and installed grackle with the cmake build system.

To help support the usage of Grackle in a wide variety of scenarios, CMake-driven installations of Grackle come with a file called **grackle.pc**\ .
For the uninitiated, this file encodes a variety of metadata, including Grackle's version number, the compilation requirements, and the linking requirements in a standardized file format.
This is the most commonly used format for specifying linking requirements on posix operating systems (including Linux, macOS, the BSDs, etc.).
The format is understood by `pkg-config <https://www.freedesktop.org/wiki/Software/pkg-config/>`__  (or an alternative implementation called `pkgconf <https://github.com/pkgconf/pkgconf>`__ ), and pkg-config comes preinstalled at many computing facilities.
This file format is recognized by most popular build systems like autotools, Meson, or even CMake (if using CMake, you should prefer the methods described :ref:`here <cmake_grackle_linking>` or :ref:`here <Embed_Grackle_in_Sim_Build>`).

If your application's build system consists of bare Makefiles, you can employ this file by invoking the ``pkg-config`` directly.
The basic usage is extremely simple:

  * ``pkg-config --cflags grackle`` provides compiler flags (namely the ``-I`` flag)

  * ``pkg-config --libs grackle`` provides linker flags (namely the ``-L`` and ``-l`` flags)

If Grackle isn't installed in a standard system installation directly, you or the end-user needs to set the ``PKG_CONFIG_PATH`` variable to tell ``pkg-config`` where to find **grackle.pc** (if Grackle is a shared library, the relevant runtime-challenges LINK still need to be addressed).

To promote a seamless user experience, the contents of **grackle.pc** are customized based on whether Grackle is installed as a shared library or as a static library.
This is the ONLY *officially* supported way to consume grackle as a static library in a non-CMake build.

The following snippet shows a sample Makefile for compiling a sample application while using Grackle.

.. code-block:: makefile

   # if Grackle is installed in an atypical location:
   # -> it is the caller's responsibility to appropriately adjust the 
   #    PKG_CONFIG_PATH environment variable so that pkg-config can find
   #    grackle.pc
   # -> it is also the the caller's responsibility to setup LD_LIBRARY_PATH
   #    appropriately if they want to use Grackle as a shared library.
   #    (Alternative extra logic can be added to add -rpath to the linker
   #    flags to accomplish the same thing)

   CFLAGS = `pkg-config --cflags grackle`
   LDFLAGS = `pkg-config --libs grackle`

   # flags unrelated to Grackle
   UNAME := $(shell uname)
   ifneq ($(UNAME), Darwin)
     OTHER_LDFLAGS=-lm
   endif

   c_example:
   	$(CC) $(CFLAGS) -c c_example.c -o c_example.o
   	$(CC) $(LDFLAGS) $(OTHER_LDFLAGS) c_example.o -o c_example

pkg-config also provides additional functionality, like querying version numbers, enforcing version requirements, etc.
Most of that functionality is described in `this guide <https://people.freedesktop.org/~dbn/pkg-config-guide.html>`__.
You can also query Grackle-related details, such as:

* the full version string (to determine if it's a dev-version or not) via ``pkg-config --variable=GRACKLE_VERSION_STR grackle``

* whether Grackle was compiled with double precision, via ``pkg-config --variable=GRACKLE_USE_DOUBLE grackle``

* whether Grackle was compiled with openmp, via ``pkg-config --variable=GRACKLE_USE_OPENMP grackle``

* the path to the :ref:`grdata cli tool <manage-data-files>` associated with this version of Grackle, via ``pkg-config --variable=GRACKLE_GRDATA_TOOL_PATH grackle`` (this might be useful for testing purposes)


.. warning::

   If the end-user uses CMake to create an installation that features Grackle as both a shared library and as a static library, we have included custom logic to try to ensure that the installed version of the **grackle.pc** file provides out-of-the-box support for the shared library version.
   This decision was made to follow established conventions.

   For properly configured files, ``pkg-config`` supports the ``--static`` flag as a way to theoretically allow downstream applications to switch between using shared and static libraries in these type of installations.
   Unfortunately, for a :ref:`variety of reasons <pkgconfig_rationale>` outside of our control, this **IS NOT** a reliable/portable solution; while it may work in some cases, it definitely won't give the desired result (or work at all) on several common platforms.
   We primarily provide this information for people who know what they are doing and want to programatically construct compiler flags for static linking based on a series of ``pkg-config`` queries.

.. note::

   At this time, pkg-config will **ONLY** work with a complete Grackle installation (i.e., it won't work with linking Grackle from a build directory).

   In the future, we may add support for creating a **grackle-uninstalled.pc** file to support linking against Grackle when it is in the build directory.


.. _cmake_grackle_linking:

CMake's ``find_package``
------------------------


.. note::

   This approach **ONLY** works if the end-user built and installed Grackle with the CMake build-system.

   We have also added experimental support for using this approach with a build directory.

CMake builds of Grackle install a Package Config File alongside the Grackle library that assists with importing information about an installation into your CMake project when you call the ``find_package`` command.
Here is a sample snippet showing how this works

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.16)
   project(GrackleExample LANGUAGES C Fortran)

   find_package(Grackle 3.3.1 REQUIRED)

   add_executable(example_app src/c_example.c)
   target_link_libraries(example_app Grackle::Grackle)

If Grackle is installed in a non-standard location, hints about its location can be specified with the ``Grackle_ROOT`` variable (or some other variables).

The logic has has been customized for the case when :ref:`shared and copies of Grackle are both installed <cmake_shared_and_static>` (it is inspired by behavior from hdf5).
``find_package`` will only import one of those libraries and it will import it as the ``Grackle::Grackle`` target.

* The caller can express a preference by requesting  ``shared`` or ``static`` component from `find_package <https://cmake.org/cmake/help/latest/command/find_package.html>`__.
  If the preference is listed after the ``COMPONENTS`` keyword, ``find_package`` considers the request to be a strong requirement (``find_package`` reports a failure if the requested type isn't installed).
  If the preference is listed after the ``OPTIONAL_COMPONENTS`` keyword, then the request is considered a weak preference (``find_package`` imports the non-preferred option if that is the only available choice).

* If the caller doesn't express any preference, a weak preference is inferred based on the current value of the ``BUILD_SHARED_LIBS`` variable.

We also encode extra metadata about the Grackle build and how it was configured as custom properties on the ``Grackle::Grackle`` target.
These can be accessed with the `get_target_property <https://cmake.org/cmake/help/latest/command/get_target_property.html>`__ command.
These properties include:

* ``GRACKLE_VERSION_STR`` -- stores the full version string (including any ``-dev`` suffix)
* ``GRACKLE_USE_DOUBLE`` -- stores whether Grackle was compiled with single or double precision
* ``GRACKLE_USE_OPENMP`` -- stores whether Grackle was compiled with OpenMP

Information about the :ref:`grdata cli tool <manage-data-files>` tool that is created and built alongside this version of Grackle is exposed via the ``Grackle::grcli`` executable target.
This can be useful for testing purposes.

.. _Embed_Grackle_in_Sim_Build:

Embedding Grackle into your Simulation Build
--------------------------------------------

If your simulation code is built with CMake, this is arguably the most convenient choice for your end-users.
Essentially, the idea is that you are compiling Grackle directly as part of your simulation.
In a sense you are providing automatic dependency management.
You can do this by making Grackle a git-submodule or using CMake's ``FetchContent`` machinery.

Be aware that if your simulation code doesn't use Fortran, you will need to add ``Fortran`` to the top-level ``project`` command OR call ``enable_language(Fortran)`` in your simulation's top level ``CMakeLists.txt`` file.
If you don't do this, linking errors can arise in certain scenarios. [#f1]_

Here are some basic code snippets showing the 2 approaches.
For simplicity, we assume Grackle is a required dependency:

1. This first snippet shows a case with git-submodule

   .. code-block:: cmake

      cmake_minimum_required(VERSION 3.16)
      project(GrackleExample LANGUAGES C Fortran)

      set(GRACKLE_SUBMODULE_PATH path/to/grackle/submodule)
      if (NOT EXISTS "${GRACKLE_SUBMODULE_PATH}")
        message(FATAL_ERROR "you forgot to initialize the Grackle submodule")
      endif()

      # configure your grackle build
      set(GRACKLE_USE_DOUBLE ON)
      set(GRACKLE_USE_OPENMP OFF)

      add_subdirectory("${GRACKLE_SUBMODULE_PATH}")

      add_executable(example_app src/c_example.c)
      target_link_libraries(example_app Grackle::Grackle)

2. This second snippet shows a case with ``FetchContent``

   .. code-block:: cmake

      cmake_minimum_required(VERSION 3.16)
      project(GrackleExample LANGUAGES C Fortran)

      include(FetchContent)

      # note: it's better to specify the actual commit-hash than a version
      #       tag (otherwise cmake will do a lot of extra work)
      FetchContent_Declare(Grackle
        GIT_REPOSITORY https://github.com/mabruzzo/grackle
        GIT_TAG 689be185ac55dba098309e2da9d6acdda37d1923
      )

      # configure your grackle build
      set(GRACKLE_USE_DOUBLE ON)
      set(GRACKLE_USE_OPENMP OFF)

      # download Grackle and trigger the build
      FetchContent_MakeAvailable(Grackle)

      add_executable(example_app src/c_example.c)
      target_link_libraries(example_app Grackle::Grackle)

Care has been taken while designing the CMake build-system to ensure that the ``Grackle::Grackle`` CMake target looks and acts the same regardless of whether it was produced with this strategy (embedding Grackle into your simulation code's build system) or imported via ``find_package`` (as discussed :ref:`here <cmake_grackle_linking>`).
In both cases, the target provides the same custom properties to describe information about the build.
See the :ref:`section <cmake_grackle_linking>` about ``find_package`` for more details.

Additionally, information about the :ref:`grdata cli tool <manage-data-files>` tool that is created and built alongside this version of Grackle is exposed via the ``Grackle::grcli`` executable target.

.. rubric:: Footnotes

.. [#f1] This is required by CMake.
         While we could implement some workarounds into Grackle's CMakeLists.txt files, this may not be a good idea.
         A post `has been created on the CMake forum <https://discourse.cmake.org/t/conventions-for-linking-implicit-dependencies-of-an-embedded-multi-language-static-library/11073?u=mabruzzo>`__ to solicit feedback on this topic.
         In the event that top-level project embeds both Grackle and some other CMake-project with Fortran source-code, it's best that the top-level project calls ``enable_langugage(Fortran)`` to ensure that both Grackle and the other CMake project use the same Fortran compiler.
