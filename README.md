# Grackle

[![Users' Mailing List](https://img.shields.io/badge/Users-List-lightgrey.svg)](https://groups.google.com/forum/#!forum/grackle-cooling-users)
[![CircleCI](https://circleci.com/gh/grackle-project/grackle/tree/main.svg?style=svg)](https://circleci.com/gh/grackle-project/grackle/tree/main)
[![Documentation Status](https://readthedocs.org/projects/grackle/badge/?version=latest)](https://grackle.readthedocs.io/en/latest/?badge=latest)

[Documentation](https://grackle.readthedocs.io/en/latest/) |
[Installation](https://grackle.readthedocs.io/en/latest/Installation.html) |
[Pygrackle Installation](https://grackle.readthedocs.io/en/latest/Python.html) |
[Usage Guide](https://grackle.readthedocs.io/en/latest/Interaction.html) |
[Integration Guide](https://grackle.readthedocs.io/en/latest/Integration.html) |
[Contributing](https://grackle.readthedocs.io/en/latest/Contributing.html) |
[Getting Help](https://grackle.readthedocs.io/en/latest/Help.html)

Grackle is a chemistry and radiative cooling library for astrophysical simulations and models.
The core library provides interfaces for C, C++ and Fortran simulation codes.
The project also offers the Pygrackle package to provide Python bindings.

## Features

Grackle provides functions to update chemistry species; solve radiative
cooling and update internal energy; and calculate cooling time, temperature,
pressure, and ratio of specific heats (γ).
The library offers

- two options for primordial chemistry and cooling. It can (i) evolve a non-equilibrium chemistry network  **OR** (ii) use tabulated cooling rates calculated with the photo-ionization code, [Cloudy](http://nublado.org).

- tabulated metal cooling rates calculated with [Cloudy](http://nublado.org).

- photo-heating and photo-ionization (with optional self-shielding corrections) from either the [Faucher-Giguere et al. (2009)](http://adsabs.harvard.edu/abs/2009ApJ...703.1416F) or [Haardt & Madau (2012)](http://adsabs.harvard.edu/abs/2012ApJ...746..125H) UV backgrounds.

- support for user-provided arrays of volumetric and specific heating rates.

Our [method paper](http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S) provides more information.

## Projects that provide out-of-the-box support for Grackle

Grackle is a popular tool (the [method paper](http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S) has over 275 citations) and has been used in a wide variety of calculations.
Below, we list open source projects that provide out-of-the-box support for Grackle:

[ChaNGa](https://github.com/N-BodyShop/changa),
[Cholla](https://github.com/cholla-hydro/cholla),
[Enzo](https://enzo-project.org/),
[Enzo-E](https://enzo-e.readthedocs.io/en/latest/),
[Gamer](https://github.com/gamer-project/gamer),
[Gasoline](https://github.com/N-BodyShop/gasoline),
[GIZMO](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html),
[Swift](https://github.com/SWIFTSIM/SWIFT)

> We welcome PRs to add your simulation code (or python package) to this list.

## Getting Grackle

Currently, Grackle must be built from source.
If you only need Grackle as a dependency of a simulation code and that code is built with CMake, then that code's build system might be configured to automatically fetch, build, and link Grackle into the code for you.

If you contribute to a simulation code, our [Integration Guide](https://grackle.readthedocs.io/en/latest/Integration.html) provides guidance on simplifying the process (for you and your users) of configuring your code to use Grackle.

### Building the Core Grackle Library From Source

Grackle requires a C99 compiler, a Fortran compiler, and HDF5 (1.6 or newer).
On most platforms, compilation with the CMake build system (3.16 or newer) is as simple as:

```shell
cmake -B build          # configure the build-directory
cmake --build ./build   # perform the build
```

You can run invoke the examples from within the ``build/examples`` directory.[^1]
To build Grackle as a shared lib, replace the first command with ``cmake -DBUILD_SHARED_LIBS=ON -Bbuild``.
To install Grackle, invoke ``cmake --install ./build [--prefix <prefix/path>]`` (the optional part lets you specify an install-path).

For more details **(especially if you encounter any errors),** see our comprehensive [Installation Guide](https://grackle.readthedocs.io/en/latest/Installation.html).
It provides more context for inexperienced CMake users, describes additional configuration options (relevant if you encounter issues), and describes Grackle's "classic" build-system.

### Building Pygrackle from Source

Once you have a Fortran compiler and a copy of HDF5 (1.6 or newer), simply invoke the following from the root of the Grackle repository

```shell
~/grackle $ pip install .
```

For more about installation see our [Pygrackle installation guide](https://grackle.readthedocs.io/en/latest/Python.html).

## Getting Started

To help you start using Grackle, we provide:

- a [Usage Guide](https://grackle.readthedocs.io/en/latest/Interaction.html)
- example Grackle programs written in [C, C++, and Fortran](https://github.com/grackle-project/grackle/tree/main/src/example)
- an [Integration Guide](https://grackle.readthedocs.io/en/latest/Integration.html) (for linking Grackle)
- a curated [guide](https://grackle.readthedocs.io/en/latest/Python.html#running-the-example-scripts) for the Pygrackle examples

## Contributing

Grackle is a community project!
We welcome patches, features, and bugfixes from any member of the community!
For more details, please see our [Constribution Guide](https://grackle.readthedocs.io/en/latest/Contributing.html) and our [Code of Conduct](https://grackle.readthedocs.io/en/latest/Conduct.html)


[^1]: You currently **NEED** to invoke the examples from the directory where they are located.
      We have a ["fix" in the pipeline](https://github.com/grackle-project/grackle/pull/246) to make this more flexible.
