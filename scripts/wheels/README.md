This directory holds scripts/resources used for precompiling wheels. They are all launched by cibuildwheel.

## Background

To precompile wheels, we make use of [cibuildwheel](https://cibuildwheel.pypa.io/en/stable/), which is a tool maintained by the [Python Packaging Authority](https://www.pypa.io/en/latest/). ``cibuildwheel`` is an extremely popular tool used for creating binary wheels (e.g. numpy, scipy, h5py, yt, pandas, astropy).

cibuildwheel is run by GitHub Actions, and its configuration values are stored within **pyproject.toml**.

## About the scripts

Gracklepy is similar to packages like numpy/scipy/h5py in the sense that its extension modules rely upon external shared libraries that aren't are always provided by external platforms (or external platforms may provide incompatible versions of the dependencies). Consequently, we need to retrieve/compile external dependencies and distribute precompiled-copies of these dependencies as part of the binary wheel. Currently, we redistribute libhdf5, HDF5's runtime dependencies, and gfortran's runtime libraries (the precise details vary with platform). Of course, we also need to distribute the associated licenses in the binary wheel.

To properly configure wheels to do all of this, we instruct cibuildwheel to invoke the scripts in this directory. Specifically, the scripts with the prefix "cibw_" are the ones that are directly invoked from cibuildwheel.
