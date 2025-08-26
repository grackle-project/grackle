Grackle
=======


.. image:: https://img.shields.io/badge/Users-List-lightgrey.svg
   :target: https://groups.google.com/forum/#!forum/grackle-cooling-users

.. image:: https://circleci.com/gh/grackle-project/grackle/tree/main.svg?style=svg
   :target: https://circleci.com/gh/grackle-project/grackle/tree/main

.. image:: https://results.pre-commit.ci/badge/github/grackle-project/grackle/main.svg
   :target: https://results.pre-commit.ci/latest/github/grackle-project/grackle/main

.. image:: https://readthedocs.org/projects/grackle/badge/?version=latest
   :target: https://grackle.readthedocs.io/en/latest/?badge=latest)

.. image:: https://img.shields.io/pypi/v/gracklepy?label=gracklepy%40pypi
   :target: https://pypi.org/project/gracklepy/

.. image:: https://img.shields.io/pypi/pyversions/gracklepy
   :target: https://pypi.org/project/gracklepy/

`Documentation <https://grackle.readthedocs.io/en/latest/>`__ |
`Installation <https://grackle.readthedocs.io/en/latest/Installation.html>`__ |
`Gracklepy Installation <https://grackle.readthedocs.io/en/latest/Python.html>`__ |
`Usage Guide <https://grackle.readthedocs.io/en/latest/Interaction.html>`__ |
`Integration Guide <https://grackle.readthedocs.io/en/latest/Integration.html>`__ |
`Contributing <https://grackle.readthedocs.io/en/latest/Contributing.html>`__ |
`Getting Help <https://grackle.readthedocs.io/en/latest/Help.html>`__

.. COMMENT:  README-MAIN-BODY-START-ANCHOR

Grackle is a chemistry and radiative cooling library for astrophysical simulations and models.
The core library provides interfaces for C, C++ and Fortran simulation codes.
The project also offers the Gracklepy package to provide Python bindings.

Features
--------

Grackle provides functions to update chemistry species; solve radiative
cooling and update internal energy; and calculate cooling time, temperature,
pressure, and ratio of specific heats (γ).
The library offers

- two options for primordial chemistry and cooling. It can (i) evolve a non-equilibrium chemistry network  **OR** (ii) use tabulated cooling rates calculated with the photo-ionization code, `Cloudy <http://nublado.org>`__.

- tabulated metal cooling rates calculated with `Cloudy <http://nublado.org>`__.

- photo-heating and photo-ionization (with optional self-shielding corrections) from either the `Faucher-Giguere et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...703.1416F>`__ or `Haardt & Madau (2012) <http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`__ UV backgrounds.

- support for user-provided arrays of volumetric and specific heating rates.

Our `method paper <http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`__ provides more information.

Projects that provide out-of-the-box support for Grackle
--------------------------------------------------------

Grackle is a popular tool (the `method paper <http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`__ has over 300 citations) and has been used in a wide variety of calculations.
Below, we list open source projects that provide out-of-the-box support for Grackle:

`ChaNGa <https://github.com/N-BodyShop/changa>`__,
`Cholla <https://github.com/cholla-hydro/cholla>`__,
`Enzo <https://enzo-project.org/>`__,
`Enzo-E <https://enzo-e.readthedocs.io/en/latest/>`__,
`Gamer <https://github.com/gamer-project/gamer>`__,
`Gasoline <https://github.com/N-BodyShop/gasoline>`__,
`GIZMO <http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html>`__,
`Swift <https://github.com/SWIFTSIM/SWIFT>`__

   We welcome PRs to add your simulation code (or python package) to this list.

Getting Grackle
---------------

Currently, the core Grackle library must be built from source.
If you only need Grackle as a dependency of a simulation code and that code is built with CMake, then that code's build system might be configured to automatically fetch, build, and link Grackle into the code for you (`Enzo-E <https://enzo-e.readthedocs.io/en/latest/>`__ is an example of a code configured in this manner).

If you contribute to a simulation code, our `Integration Guide <https://grackle.readthedocs.io/en/latest/Integration.html>`__ provides guidance on simplifying the process (for you and your users) of configuring your code to use Grackle.

Building the Core Grackle Library From Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grackle requires a C99 compiler, a Fortran compiler, and HDF5 (1.6 or newer).
On most platforms, compilation with the CMake build system (3.16 or newer) is as simple as:

.. code-block:: shell-session

   cmake -B build          # configure the build-directory
   cmake --build ./build   # perform the build

You can invoke the examples from within the ``build/examples`` directory (`we're adding support <https://github.com/grackle-project/grackle/pull/246>`__ to let you run them from anywhere).
To build Grackle as a shared lib, replace the first command with ``cmake -DBUILD_SHARED_LIBS=ON -Bbuild``.
To install Grackle, invoke ``cmake --install ./build [--prefix <prefix/path>]`` (the optional part lets you specify an install-path).

For more details **(especially if you encounter any errors),** see our comprehensive `Installation Guide <https://grackle.readthedocs.io/en/latest/Installation.html>`__.
It provides more context for inexperienced CMake users, describes additional configuration options (relevant if you encounter issues), and describes Grackle's "classic" build-system.

Installing Gracklepy
~~~~~~~~~~~~~~~~~~~~

The easiest way to get Gracklepy is to invoke the following command

.. code-block:: shell-session

   ~/grackle $ pip install gracklepy

This will download a pre-built version (called a wheel) of Gracklepy from PyPI. This should "just work," unless you use a highly unusual system (if it fails please let us know).

Be aware that the vast majority of Grackle calculation requires Grackle's data files.
At this time you must download manually download these files (see the installation guide); we're working on streamlining this in the future.

For more about installation (and downloading data files), see our `Gracklepy installation guide <https://grackle.readthedocs.io/en/latest/Python.html>`__ .

**NOTE:** Gracklepy was formerly known as Pygrackle.
If you previously installed Pygrackle, you should uninstall it before you install GracklePy.

Getting Started
---------------

To help you start using Grackle, we provide:

- a `Usage Guide <https://grackle.readthedocs.io/en/latest/Interaction.html>`__
- example Grackle programs written in `C, C++, and Fortran <https://github.com/grackle-project/grackle/tree/main/src/example>`__
- an `Integration Guide <https://grackle.readthedocs.io/en/latest/Integration.html>`__ (for linking Grackle)
- a curated `guide <https://grackle.readthedocs.io/en/latest/Python.html#running-the-example-scripts>`__ for the Gracklepy examples

Contributing
------------

Grackle is a community project!
We welcome patches, features, and bugfixes from any member of the community!
For more details, please see our `Constribution Guide <https://grackle.readthedocs.io/en/latest/Contributing.html>`__ and our `Code of Conduct <https://grackle.readthedocs.io/en/latest/Conduct.html>`__

Citing Grackle
--------------

If you use Grackle please cite it.
More instructions are provided `here <https://grackle.readthedocs.io/en/latest/Citing.html>`__.
