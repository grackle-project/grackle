# Grackle

[![Users' Mailing List](https://img.shields.io/badge/Users-List-lightgrey.svg)](https://groups.google.com/forum/#!forum/grackle-cooling-users)
[![Build Status](https://travis-ci.org/grackle-project/grackle.svg?branch=master)](https://travis-ci.org/grackle-project/grackle)
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

## Resources

- documentation: https://grackle.readthedocs.io/

- source code repository: https://github.com/grackle-project/grackle

- method paper: [Smith et al. (2017)](http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S)
