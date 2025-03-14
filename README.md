> [!IMPORTANT]  
> This is a special version of the [brittonsmith:gen2024](https://github.com/brittonsmith/grackle/tree/gen2024) branch (i.e. the branch of changes proposed for merging in the grackle-project/grackle#177 Pull Request).
>
> This branch includes additional changes that are needed to simplify the transcription process to C++. There are pending PRs to merge all of these changes into the gen2024 branch ([see this list of PRs](https://github.com/brittonsmith/grackle/pulls/mabruzzo))


# Grackle

[![Users' Mailing List](https://img.shields.io/badge/Users-List-lightgrey.svg)](https://groups.google.com/forum/#!forum/grackle-cooling-users)
[![CircleCI](https://circleci.com/gh/grackle-project/grackle/tree/main.svg?style=svg)](https://circleci.com/gh/grackle-project/grackle/tree/main)
[![Documentation Status](https://readthedocs.org/projects/grackle/badge/?version=latest)](https://grackle.readthedocs.io/en/latest/?badge=latest)

Grackle is a chemistry and radiative cooling library for astrophysical
simulations and models.  Grackle has interfaces for C, C++, Fortran, and
Python codes and provides:

- two options for primordial chemistry and cooling:

   1. non-equilibrium primordial chemistry network for atomic H, D, and He
   as well as H<sub>2</sub> and HD, including H<sub>2</sub> formation on dust grains.

   2. tabulated H and He cooling rates calculated with the photo-ionization
      code, [Cloudy](http://nublado.org).

- tabulated metal cooling rates calculated with [Cloudy](http://nublado.org).

- photo-heating and photo-ionization from two UV backgrounds with optional
  self-shielding corrections:

   1. [Faucher-Giguere et al. (2009)](http://adsabs.harvard.edu/abs/2009ApJ...703.1416F).

   2. [Haardt & Madau (2012)](http://adsabs.harvard.edu/abs/2012ApJ...746..125H).

- support for user-provided arrays of volumetric and specific heating rates.

The Grackle provides functions to update chemistry species; solve radiative
cooling and update internal energy; and calculate cooling time, temperature,
pressure, and ratio of specific heats (gamma).

For more information on features, installation, and integration with simulation
codes and models, see our [online documentation](https://grackle.readthedocs.io/).

## Software that uses Grackle

A (non-exhaustive) list of software that provides out-of-the-box support for using Grackle includes:

- [ChaNGa](https://github.com/N-BodyShop/changa)

- [Cholla](https://github.com/cholla-hydro/cholla)

- [Enzo](https://enzo-project.org/)

- [Enzo-E](https://enzo-e.readthedocs.io/en/latest/)

- [Gamer](https://github.com/gamer-project/gamer)

- [Gasoline](https://github.com/N-BodyShop/gasoline)

- [GIZMO](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html)

- [Swift](https://github.com/SWIFTSIM/SWIFT)

We welcome PRs to add your simulation code to this list. We also welcome the inclusion of python modules that depend on Pygrackle.

## Resources

- documentation: https://grackle.readthedocs.io/

- source code repository: https://github.com/grackle-project/grackle

- method paper: [Smith et al. (2017)](http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S)
