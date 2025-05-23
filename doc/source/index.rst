.. grackle documentation master file, created by
   sphinx-quickstart on Sun Mar 24 09:23:44 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to grackle's documentation!
===================================

Grackle is a chemistry and radiative cooling library for astrophysical simulations and models.
The core library provides interfaces for C, C++ and Fortran simulation codes.
The project also offers the Pygrackle package to provide Python bindings.

Features
--------
Grackle provides functions to update chemistry species; solve radiative
cooling and update internal energy; and calculate cooling time, temperature,
pressure, and ratio of specific heats (γ).
The library offers

- two options for primordial chemistry and cooling. It can (i) evolve a non-equilibrium chemistry network  **OR** (ii) use tabulated cooling rates calculated with the photo-ionization code, `Cloudy <http://nublado.org>`_.

- tabulated metal cooling rates calculated with `Cloudy <http://nublado.org>`_.

- photo-heating and photo-ionization (with optional self-shielding corrections) from either the `Faucher-Giguere et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...703.1416F>`_ or `Haardt & Madau (2012) <http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`_ UV backgrounds.

- support for user-provided arrays of volumetric and specific heating rates.

Our `method paper <http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`__ provides more information.

Getting Grackle
---------------

Currently, Grackle must be built from source.
If you only need Grackle as a dependency of a simulation code and that code is built with CMake, then that code's build system might be configured to automatically fetch, build, and link Grackle into the code for you (`Enzo-E <https://enzo-e.readthedocs.io/en/latest/>`__ is an example of a code configured in this manner).

If you contribute to a simulation code, our `Integration Guide <https://grackle.readthedocs.io/en/latest/Integration.html>`__ provides guidance on simplifying the process (for you and your users) of configuring your code to use Grackle.

Building the Core Grackle Library From Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grackle requires a C99 compiler, a Fortran compiler, and HDF5 (1.6 or newer).
On most platforms, compilation with the CMake build system (3.16 or newer) is as simple as:

.. code-block:: shell-session

   cmake -B build          # configure the build-directory
   cmake --build ./build   # perform the build

You can run invoke the examples from within the ``build/examples`` directory.\ [#f1]_
To build Grackle as a shared lib, replace the first command with ``cmake -DBUILD_SHARED_LIBS=ON -Bbuild``.
To install Grackle, invoke ``cmake --install ./build [--prefix <prefix/path>]`` (the optional part lets you specify an install-path).

For more details **(especially if you encounter any errors),** see our comprehensive `Installation Guide <https://grackle.readthedocs.io/en/latest/Installation.html>`__.
It provides more context for inexperienced CMake users, describes additional configuration options (relevant if you encounter issues), and describes Grackle's "classic" build-system.

Building Pygrackle from Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have a Fortran compiler and a copy of HDF5 (1.6 or newer), simply invoke the following from the root of the Grackle repository

.. code-block:: shell-session

   ~/grackle $ pip install .

For more about installation see our `Pygrackle installation guide <https://grackle.readthedocs.io/en/latest/Python.html>`__.

Getting Started
---------------

To help you start using Grackle, we provide:

- a `Usage Guide <https://grackle.readthedocs.io/en/latest/Interaction.html>`__
- example Grackle programs written in `C, C++, and Fortran <https://github.com/grackle-project/grackle/tree/main/src/example>`__
- an `Integration Guide <https://grackle.readthedocs.io/en/latest/Integration.html>`__ (for linking Grackle)
- a curated `guide <https://grackle.readthedocs.io/en/latest/Python.html#running-the-example-scripts>`__ for the Pygrackle examples

.. toctree::
   :hidden:
   :maxdepth: 2

   Installation.rst
   Testing.rst
   Tutorial.rst
   Integration.rst
   Interaction.rst
   Parameters.rst
   RateFunctions.rst
   Reference.rst
   Versioning.rst
   Python.rst
   Conduct.rst
   Contributing.rst
   Help.rst

.. include:: Help.rst

Contributing
------------

Development of Grackle happens in the open on GitHub `here
<https://github.com/grackle-project/grackle>`__.  We welcome new
contributors. Please, see the :ref:`conduct`.  For a guide to developing
Grackle, see :ref:`contributing-code`.

Citing grackle
--------------

The Grackle method paper can be found
`here <http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`__.

The Grackle library was born out of the chemistry and cooling routines of the 
`Enzo <http://enzo-project.org/>`__ simulation code.  As such, all of those who 
have contributed to Enzo development, and especially to the chemistry and 
cooling, have contributed to the Grackle.

If you used the Grackle library in your work, please cite it as "the Grackle 
chemistry and cooling library (`Smith et al. 2017
<http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`_)."  Also, please add a
footnote to `https://grackle.readthedocs.io/ <https://grackle.readthedocs.io/>`_.

.. rubric:: Footnotes

.. [#f1] You currently **NEED** to invoke the examples from the directory where they are located.
         We have a `"fix" in the pipeline <https://github.com/grackle-project/grackle/pull/246>`__ to make this more flexible.
