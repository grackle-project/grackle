.. _python:

Pygrackle: Running Grackle in Python
====================================

Grackle comes with a Python interface, called Pygrackle, which provides
access to all of Grackle's functionality.

To install Pygackle, you'll need to make sure that HDF5 and a fortran compiler are installed (for building the Grackle library itself).

Pygrackle's runtime-dependencies are:

- `h5py <https://www.h5py.org/>`__
- `matplotlib <https://matplotlib.org/>`__
- `NumPy <https://www.numpy.org/>`__
- `yt <https://yt-project.org/>`__

The above dependencies are automatically installed with pip (alternatively you can follow instructions for installing yt).

If you want to run the pygrackle test-suite, you'll need the ``packaging`` and ``py.test`` packages. If pip is up to date (25.1 or newer), you can simply invoke ``pip install --group dev`` from the root of the repository.

.. _install-pygrackle:

Installing Pygrackle
--------------------

Currently, the only way to get Pygrackle is to build it from source.

There are 3 ways to build Pygrackle:

1. As a standalone, self-contained module.
   The build command creates a fresh build of the Grackle library and packages it with the Pygrackle module.
   **(This is the recommended approach)**

2. As a module that links to an external copy of Grackle that was compiled with the :ref:`Classic build system <classic_build>`.
   (This is consistent with the legacy approach for building Pygrackle).

3. As a module that links to an external copy of Grackle that was created with the :ref:`CMake build system <cmake_build>`.

Currently, Pygrackle should be used with Grackle builds where OpenMP is disabled.

.. warning::

   We strongly encourage you to use the first approach so that your Pygrackle installation is independent of other Grackle installations on your machine.

   The latter 2 approaches are primarily intended for testing-purposes.
   If you use the latter 2 approaches, it's your responsibility to ensure that you delete the old version of Pygrackle and reinstall it whenever the external Grackle library is updated.
   If you forget, Pygrackle may still work, but it's more likely to produce a segmentation fault or (even worse!) silently give incorrect results.

1. Build Pygrackle as a standalone module (recommended)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

To install Pygrackle, you just need to invoke the following from the root directory.

.. code-block:: shell-session

    ~/grackle $ pip install -v .

Under this approach, pygrackle's python build backend, `scikit-build-core <https://scikit-build-core.readthedocs.io/en/latest/>`__, automatically builds the core Grackle library (:ref:`in a CMake build <cmake_build>`), packages the resulting library as part of the pygrackle package, and cleans up from the build.
If a new enough version of CMake cannot be found, the above command will automatically download CMake (it will only be accessible during the build).

While the above command should "just work," you have a few options for customizing the build:

1. In most cases, you probably want to use the standard environment variables directly understood by CMake:\ [#f1]_

   - you can specifiy your choice of C and Fortran compilers that are used for this by manipulating the :envvar:`!CC` and :envvar:`!FC` environment variables.
   - if you must pass extra compiler flags to all invocations of the C or Fortran compiler (this *shouldn't* really be necessary), you can use the :envvar:`!CFLAGS` or :envvar:`!FFLAGS` environment variable.
   - to provide a hint about the location of the hdf5 library you can use :envvar:`!HDF5_ROOT` or :envvar:`!HDF5_DIR` (these behave similarly to the CMake variables of the same name).

2. Alternatively, you can specify the values of CMake variables.

   - Recall that when you are directly using the CMake build system to build a project, you configure build-properties by defining CMake variables.
     :ref:`As we explain elsewhere, <how_to_configure>` this is commonly accomplished by listing arguments of the form ``-D<variable>=<value>`` when we call CMake on the command line.
   - Things are a little different in the context of building pygrackle, since we aren't  invoking CMake, directly.
     Instead we need to instruct the `scikit-build-core <https://scikit-build-core.readthedocs.io/en/latest/>`__ python build-backend to forward arguments onto CMake.
   - There are effectively 2 ways of doing this: (i) specify this information as extra command line arguments when invoking ``pip install`` OR (ii) we can specify it through an environment variable understood by scikit-build-core

   For the sake of concreteness, let's imagine that you want to assign the CMake variables ``CMAKE_C_COMPILER`` and ``HDF5_ROOT`` values of ``gcc-14`` and ``/path/to/hdf5``.
   The following code-snippets illustrate how to do this:

   .. tabs::

      .. group-tab:: pip argument

         If using a version of pip from before 23.1, you need to replace ``-C`` with ``--config-settings=`` in the following snippet.

         .. code-block:: shell-session

            ~/grackle $ pip install -v . -Ccmake.args=-DCMAKE_C_COMPILER=gcc-14;-DHDF5_ROOT=/path/to/hdf5

      .. group-tab:: Environment Variable

         You can use the use the :envvar:`!SKBUILD_CMAKE_ARGS` env variable, where each argument is separated by semicolons, or :envvar:`!CMAKE_ARGS`, where each argument is separated by a space.

         .. tabs::

            .. code-tab:: shell-session :envvar:`!SKBUILD_CMAKE_ARGS`

               ~/grackle $ export SKBUILD_CMAKE_ARGS="-DCMAKE_C_COMPILER=gcc-14;-DHDF5_ROOT=/path/to/hdf5"
               ~/grackle $ pip install -v .


            .. code-tab:: shell-session :envvar:`!CMAKE_ARGS`

               ~/grackle $ export CMAKE_ARGS="-DCMAKE_C_COMPILER=gcc-14 -DHDF5_ROOT=/path/to/hdf5"
               ~/grackle $ pip install -v .


.. tip::

   If you are looking to modify a standard CMake option, you should generally check scikit-build-core's `documentation  <https://scikit-build-core.readthedocs.io/en/latest/>`__; there are some special cases.

   Consider the ``CMAKE_BUILD_TYPE`` variable, which controls :ref:`optimiziation flags and the presence of debugger symbols <how_to_configure>`.
   Rather than directly modifying this variable, you should modify scikit-build-core's ``cmake.build-type`` variable.
   If you wanted to set it to ``Debug``, you might do one of the following options:

   .. tabs::

      .. group-tab:: pip argument

         If using a version of pip from before 23.1, you need to replace ``-C`` with ``--config-settings=`` in the following snippet.

         .. code-block:: shell-session

            ~/grackle $ pip install -v . -Ccmake.build-type="Debug"

      .. group-tab:: Environment Variable

         You can would store ``"Debug"`` within the :envvar:`!SKBUILD_CMAKE_BUILD_TYPE` env variable


If you encounter any compilation problems, you can also link Pygrackle against a version of the Grackle library that you already built.

(In the event that you are writing an external python package that depends on directly linking to the underlying Grackle library, be aware that the underlying organization of files in the resulting package may change.
We have no plans to support this scenario.)

2. Link to external Grackle library (built with Classic build system)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*Prerequisite:* This scenario assumes that you have already built Grackle with the :ref:`Class build system <classic_build>`.

To build Pygrackle, we use the ``PYGRACKLE_LEGACY_LINK`` environment variable to indicate that we want to that external build.
Specifically, that variable must be configured as ``PYGRACKLE_LEGACY_LINK=classic``.

.. code-block:: shell-session

    ~/grackle $ PYGRACKLE_LEGACY_LINK=classic pip install .

.. note::

   We explicitly try to maintain the legacy behavior of the older setuptools-based python build-system.
   This means that we use the copy of the Grackle shared library from the build directory during linking (i.e. Pygrackle will happily build even if Grackle isn't fully installed).
   
   We then **ASSUME** that a copy of the Grackle shared library will be in a location known to the system, when you try to run Pygrackle.
   This could be a standard system location for libraries (on some systems you may need to invoke ``ldconfig`` after installation).
   This could also be a location specified by the relevant variable; ``LD_LIBRARY_PATH`` if you're on Linux (or most unix-like systems) or ``DYLD_LIBRARY_PATH`` (if you're on macOS)

3. Link to external Grackle library (built with the CMake build system)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*Prerequisite:*  This scenario assumes that you have already built (and possibly installed) Grackle with the :ref:`CMake build system <cmake_build>`.
Specifically, that cmake build must have compiled Grackle as a shared library (the primary way to ensure this happens is by passing the ``-DBUILD_SHARED_LIBS=ON`` flag when using ``cmake`` to configure the build).

To build Pygrackle in this way, you must initialize either the ``Grackle_DIR`` environment variable or the ``Grackle_ROOT`` environment variable with the relevant path for your prebuilt Grackle library.
This path can either point to cmake build directory (where Grackle is built) OR an installation directory.

We illustrates how to install Pygrackle under this approach down below.
For the sake of example, we assume that we previously used ``cmake`` to build (and compile) Grackle as a shared library in a build directory called **~/grackle/build**.

.. tabs::

   .. tab:: Default Case (libgrackle won't move after building)

      The default command to build Pygrackle against a CMake-built is shown below.
      **By default, this approach assumes that the Grackle shared library will never move.**
      This means that issues will occur if you delete or move the Grackle library.
      (This is a necessary assumption in order to support build directories).

      .. code-block:: shell-session

         ~/grackle $ Grackle_DIR=${PWD}/build pip install .

   .. tab:: Legacy Linking

      It's also possible to achieve linking behavior more similar to the case where we build Pygrackle against an external Grackle library that was built with the classic build system (this is consistent with the behavior implemented by Pygrackle's former ``setuptools`` build system).
      Under this scenario, no relationship is assumed between the path to the Grackle shared library that is used while building Pygrackle and the path that is used while running Pygrackle.
      Instead, we assume that the Grackle shared library will be at an arbitrary location known to the system at runtime (e.g. either it's in a standard location that the OS knows to check or you use ``LD_LIBRARY_PATH``/``DYLD_LIBRARY_PATH``.

      To easily invoke this linking behavior, you can either pass an additional argument to ``pip`` or define an environment variable.

       .. tabs::

          .. code-tab:: shell-session pip

             ~/grackle $ Grackle_DIR=${PWD}/build \
             > pip install . --config-settings=cmake.define.CMAKE_SKIP_INSTALL_RPATH=TRUE"

          .. code-tab:: shell-session Environment

             ~/grackle $ export Grackle_DIR=${PWD}/build
             ~/grackle $ export SKBUILD_CMAKE_DEFINE="CMAKE_SKIP_INSTALL_RPATH=TRUE"
             ~/grackle $ pip install --user .

Testing Your Installation
-------------------------

To make sure everything is installed properly, you can try invoking pygrackle from the command line:

.. code-block:: shell-session

   $ python -c "import pygrackle"

If this command executes without raising any errors, then you have successfully installed Pygrackle.

Installing DataFiles
++++++++++++++++++++

To install the datafiles in a location usable for automatic usage in the Pygrackle examples (and tests) we recommend invoking the following command (from any directory):

.. code-block:: shell-session

   $ python -m pygrackle fetch

:ref:`This section <manage-data-files>` for more details about customizing the the location where data is stored and about managing datafiles in general.

.. _pygrackle-dev:

Installing Pygrackle Development Requirements
+++++++++++++++++++++++++++++++++++++++++++++

There are a handful of additional packages required purely for developing
Grackle. For example, these will enable :ref:`testing` and building
the documentation locally. These dependencies are specified as dependency
groups, which can be installed with pip (v25.1).
To install all of these dependencies, you can invoke

.. code-block:: shell-session

   ~/grackle $ pip install --group dev

The above command will install the dependencies independently of Pygrackle.
To install these dependencies at the same time as Pygrackle, you can replace last line of the :ref:`pygrackle installation instructions <install-pygrackle>` with:

.. code-block:: shell-session

   ~/grackle $ pip install --group=dev -e .

The above snippet, includes the optional ``-e`` flag to perform an editable-install, which is necessary to run most tests.

.. tip::

   The high level interface of the `uv python package manager <https://docs.astral.sh/uv/>`__ automatically installs the "dev" dependency-group when you install Pygrackle from source.

Running the Example Scripts
---------------------------

A number of example scripts are available in the :source:`src/python/examples`
directory.  These scripts provide examples of ways that Grackle can be
used in simplified models, such as solving the temperature evolution of
a parcel of gas at constant density or in a free-fall model.  Each example
will produce a figure as well as a dataset that can be loaded and analyzed
with `yt <http://yt-project.org/>`__.

All of the example scripts discussed below use the following line to
make a guess at where the Grackle input files are located.

.. caution::

   This snippet is **NOT** part of the public API.
   It is a short-term solution that is being used until functionality proposed by :ghpr:`237` can be reviewed.

.. code-block:: python

   from pygrackle.utilities.data_path import grackle_data_dir

Cooling Rate Figure Example
+++++++++++++++++++++++++++

This sets up a one-dimensional grid at a constant density with 
logarithmically spaced temperatures from 10 K to 10\ :sup:`9` K.  Radiative cooling 
is disabled and the chemistry solver is iterated until the species fractions 
have converged.  The cooling time is then calculated and used to compute the cooling 
rate.

.. code-block:: shell-session

   ~/grackle/src/python/examples $ python cooling_rate.py

.. image:: _images/cooling_rate.png
   :width: 500

After the script runs, and hdf5 file will be created with a similar name.  This
can be loaded in with yt.

.. code-block:: python

   >>> import yt
   >>> ds = yt.load("cooling_rate.h5")
   >>> print ds.data["temperature"]
   [  1.00000000e+01   1.09698580e+01   1.20337784e+01   1.32008840e+01, ...,
      7.57525026e+08   8.30994195e+08   9.11588830e+08   1.00000000e+09] K
   >>> print ds.data["cooling_rate"]
   [  1.09233398e-25   1.08692516e-25   1.08117583e-25   1.07505345e-25, ...,
      3.77902570e-23   3.94523273e-23   4.12003667e-23   4.30376998e-23] cm**3*erg/s


Cooling Cell Example
++++++++++++++++++++

This sets up a single grid cell with an initial density and temperature and solves 
the chemistry and cooling for a given amount of time.  The resulting dataset gives
the values of the densities, temperatures, and mean molecular weights for all times.

.. code-block:: shell-session

   ~/grackle/src/python/examples $ python cooling_cell.py

.. image:: _images/cooling_cell.png
   :width: 500

.. code-block:: python

   >>> import yt
   >>> ds = yt.load("cooling_cell.h5")
   >>> print ds.data["time"].to("Myr")
   YTArray([  0.00000000e+00,   6.74660169e-02,   1.34932034e-01, ...,
            9.98497051e+01,   9.99171711e+01,   9.99846371e+01]) Myr
   >>> print ds.data["temperature"]
   YTArray([ 990014.56406726,  980007.32720091,  969992.99066987, ...,
             9263.81515866,    9263.81515824,    9263.81515865]) K


Free-Fall Collapse Example
++++++++++++++++++++++++++

This sets up a single grid cell with an initial number density of 1 cm\ :sup:`-3`.  
The density increases with time following a free-fall collapse model.  As the density 
increases, thermal energy is added to model heating via adiabatic compression.
This can be useful for testing chemistry networks over a large range in density.

.. code-block:: shell-session

   ~/grackle/src/python/examples $ python freefall.py

.. image:: _images/freefall.png
   :width: 500

The resulting dataset can be analyzed similarly as above.

.. code-block:: python

   >>> import yt
   >>> ds = yt.load("freefall.h5")
   >>> print ds.data["time"].to("Myr")
   [   0.            0.45900816    0.91572127 ...,  219.90360841  219.90360855
     219.9036087 ] Myr
   >>> print ds.data["density"]
   [  1.67373522e-25   1.69059895e-25   1.70763258e-25 ...,   1.65068531e-12
      1.66121253e-12   1.67178981e-12] g/cm**3
   >>> print ds.data["temperature"]
   [   99.94958248   100.61345564   101.28160228 ...,  1728.89321898
     1729.32604568  1729.75744287] K

Using Grackle with yt
+++++++++++++++++++++

This example illustrates how Grackle functionality can be called using
simulation datasets loaded with `yt <https://yt-project.org/>`__ as
input. Note, below we invoke Python with the ``-i`` flag to keep the
interpreter running. The second block is assumed to happen within the
same session.

.. code-block:: shell-session

   ~/grackle/src/python/examples $ python -i yt_grackle.py

.. code-block:: python

   >>> print (sp['gas', 'grackle_cooling_time'].to('Myr'))
   [-5.33399975 -5.68132287 -6.04043746 ... -0.44279721 -0.37466095
    -0.19981158] Myr
   >>> print (sp['gas', 'grackle_temperature'])
   [12937.90890302 12953.99126155 13234.96820101 ... 11824.51319307
    11588.16161462 10173.0168747 ] K

Through ``pygrackle``, the following ``yt`` fields are defined:

- ``('gas', 'grackle_cooling_time')``
- ``('gas', 'grackle_gamma')``
- ``('gas', 'grackle_molecular_weight')``
- ``('gas', 'grackle_pressure')``
- ``('gas', 'grackle_temperature')``
- ``('gas', 'grackle_dust_temperature')``

These fields are created after calling the ``add_grackle_fields`` function.
This function will initialize Grackle with settings from parameters in the
loaded dataset. Optionally, parameters can be specified manually to override.


.. rubric:: Footnotes

.. [#f1] When you want to overwrite a configuration-default in a CMake build, we generally encourage use of CMake variables that are counterparts to environment variables (e.g. prefer passing ``-DCMAKE_C_COMPILER=<blah>`` on the command-line to exporting ``CC=<blah>``), since the former has precedence and it provides more fine-grained control (i.e. some options are only influenced by CMake variables).
         However, in this case of driving a python build (that internally creates a CMake build, installs the products, and cleans up the build directory), the use of environment variables is somewhat more common, and you usually don't need as much fine-grain control.
         Furthermore, the names of the environment variables follow fairly standard Unix conventions.
