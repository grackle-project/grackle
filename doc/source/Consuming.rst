Integrating Grackle into your Application
=========================================

**If you are simply installing Grackle because it is a dependency of a simulation-code, you can probably disregard this section.**

This section assumes that the reader wants to integrate Grackle into their simulation code.
There are broadly 3 mainstream approaches (each with a couple variations) that your simulation code could take to support Grackle:

1. :ref:`Link_Against_Grackle_Installation` -- this is the most conventional approach

2. :ref:`Link_Against_Grackle_Build` -- this is not currently supported

3. :ref:`Embed_Grackle_in_Sim_Build` -- this is most convenient for end-users (it amounts to automatic dependency management), but it requires that your simulation code is built with CMake.

.. note::

   If you want to support end-users that are installing Grackle with the classic build-system, there is :ref:`only 1 supported way <manual_grackle_linking>`
   Unless explicitly noted otherwise, the rest of this section assumes that your code will **ONLY** consume Grackle when it is built with CMake.

.. note::

   If we ignored ease-of-linking, you should almost always prefer to consume Grackle as a static library because the resulting code will be slightly faster.
   Because of Grackle's API (and the limited number of Grackle-developers), you should always recompile your application when the Grackle library is modified, even if you are consuming Grackle as a shared library (eliminating a strength of shared libraries).
   
   With that said, we will point out a particular case where Grackle is **MUCH** easier to consume if you use it as a shared library. 


.. _Link_Against_Grackle_Installation:

Link Against Grackle Installation
---------------------------------

Essentially there are 3 flavors to this approach.

.. _manual_grackle_linking:

Manually Specifying Linking Flags (when using Grackle as a shared library)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This is the historic approach that we have always supported.
If you employ this approach, you should inform your end-users that they should employ build-and-install Grackle via the classic build system or use the CMake build-system to build (and install) Grackle as a shared library.

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
   OTHER_LDFLAGS=-lm

   c_example:
           $(CC) $(CFLAGS) c_example.c -o c_example.o
           $(CC) $(LDFLAGS) $(OTHER_LDFLAGS) c_example.o -o c_example

This scenario has simple linker flags because we are linking to the shared library form of Grackle.
If you used Grackle in its static library form, you would also need to manually link to all of Grackle's dependencies (or at least the ones that don't overlap with your simulations other dependencies.
Some of these dependencies are implicit and depend on the precise choice of compiler and whether openmp is used.

.. warning::

   If a Grackle-installation using the Classic Build system made use of an incomplete machine-file, the Grackle shared-library may not be properly linked to all required runtime dependencies.
   Consequently, the simulation code may need to link against extra dependencies
   This commonly happens with the OpenMP runtime libraries.

   The CMake-builds are much more robust against these kinds of errors.

Using pkg-config
++++++++++++++++

.. note::

   `GH-#204 <https://github.com/grackle-project/grackle/pull/204>`__ will add support for this approach for **any** installation of a cmake-build of Grackle (whether you compiled Grackle as a static or a shared library).

CMake's ``find_package`` (in Config mode)
+++++++++++++++++++++++++++++++++++++++++

.. note::

   `GH-#204 <https://github.com/grackle-project/grackle/pull/204>`__ will add support for this approach for **any** installation of a cmake-build of Grackle (whether you compiled Grackle as a static or a shared library).

   The rest of this subsection has been written as though that PR is already merged.

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
  If the preference is listed after the ``OPTIONAL_COMPONENTS`` keyword, then the request is considered a weak preference (``find_package`` import the non-preferred option if that is the only available choice).

* If the caller doesn't express any preference, a weak preference is inferred based on the current value of the ``BUILD_SHARED_LIBS`` variable.

We also encode extra-metadata about the Grackle build and how it was configured as custom properties on the ``Grackle::Grackle`` target.
These can be accessed with the `get_target_property <https://cmake.org/cmake/help/latest/command/get_target_property.html>`__ command.
These properties include:

* ``GRACKLE_VERSION_STR`` -- stores the full version string (including any ``-dev`` suffix
* ``GRACKLE_USE_DOUBLE`` -- stores whether Grackle was compiled with single or double precision
* ``GRACKLE_USE_OPENMP`` -- stores whether Grackle was compiled with OpenMP

.. _Link_Against_Grackle_Build:

Link Against a Grackle Build-Directory
--------------------------------------

This is **NOT** currently supported.

.. warning::

   We will add support for using cmake's ``find_package`` to link against the contents of a build-directory in the near future.
   We may also add support for using pkg-config for the same purpose.

   We don't currently plan to support manual linking to libraries in the build directory.
   If this is something you want to be able to do, please let the developers know.
   Be advised, the organization and precise contents of the build-directory **will** change in the short-term (e.g. some "hacky," temporary choices were made to get tests running that we intend to more properly address).

.. _Embed_Grackle_in_Sim_Build:

Embed Grackle into your Simulation Build
----------------------------------------

If your simulation-code is built with CMake, this is arguably the most convenient choice for your end-users.
Essentially, the idea is that you are compiling Grackle directly as part of your simulation.
In a sense you are providing automatic dependency management.
You can do this by making Grackle a git-submodule or using CMake's ``FetchContent`` machinery.

Be aware that if your simulation code doesn't use Fortran, you will need to add ``Fortran`` to the top-level ``project`` command OR call ``enable_language(Fortran)`` in your simulation's top level ``CMakeLists.txt`` file.
If you don't do this, linking errors can arise in certain scenarios. [#f1]_

Here are some basic code-snippets showing the 2 approaches.
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

Care has been taken while designing the CMake build-system to ensure that the CMake target produced looks and acts the same regardless of whether it was imported via ``find_package`` or produced by embedding Grackle into your simulation code.
In both cases, the target provides the same custom properties to describe information about the build.


.. rubric:: Footnotes

.. [#f1] This is required by CMake.
         While we could implement some workarounds into Grackle's CMakeLists.txt files, they all involve some assumptions.
         In the event that top-level project depends embeds both Grackle and some other CMake-project with Fortran source-code, it's best that the top-level project calls ``enable_langugage(Fortran)`` to ensure that the both Grackle and the other CMake-project use the same Fortran compiler.
