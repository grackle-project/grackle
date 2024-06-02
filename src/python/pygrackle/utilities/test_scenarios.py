########################################################################
#
# utilities for setting up interesting test-problems (this is expected
# to be used for unit-testing)
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################
import numpy as np
from .physical_constants import mass_hydrogen_cgs, boltzmann_constant_cgs

def grid_initial_values(primordial_chem_flag, metal_cooling_flag, code_units,
                        target_T_cgs_arr, target_Hmassdensity_cgs_arr,
                        target_metallicity_arr = None,
                        nominal_initial_mmw = 0.6,
                        nominal_initial_gamma = 5/3.0,
                        nominal_HMassFrac = 0.76,
                        nominal_SolarMetalFractionByMass = 0.0125,
                        nominal_DeuteriumToHydrogenMassRatio = 6.8e-5):
    """
    Create a dict holding arrays of various fields used by Grackle calculations

    A major goal of this function is to be able to initialize a grid without
    requiring a fully initialized grackle instance. It's fine if we also want to
    add support for the case where grackle is initialized, but I think we should
    maintain support for this case

    Parameters
    ----------
    primordial_chem_flag : int
        Flag that will be used in the intended grackle-calculation that controls
        which primordial chemistry network is used
    metal_cooling_flag : int
        Flag that will be used in the intended grackle-calculation to enable
        metal cooling using the Cloudy tables.
    code_units : dict
        Specify the basic code_units. The only keys that are used are
        'length_units', 'density_units', 'time_units'
    target_T_cgs_arr : array_like
        Scalar or 1D array of initial temperature values with units of Kelvin.
        These temperatures are targetted, but may not exactly be hit.
    target_Hmassdensity_cgs_arr : array_like
        Scalar or 1D array of initial mass densities of Hydrogen (whether it's
        in ionized, atomic, or molecular form
    target_metallicity_arr : array_like, Optional
        Specifies the target metallicities. When unspecified it defaults to
        zero for no metal cooling and to 1 (solar metallicity) in other cases.

    Notes
    -----
    We intentionally try to avoid using internal grackle functionality in
    this particular function. We use the term, nominal, to reflect the fact
    that the values may not necessarily reflect the actual values that are
    being used
    """

    # at the moment, we assume solar metallicity. In the future, it might be
    # nice to support varying the metallicity.

    # we make a lot of simplifying assumptions just to get a decent result. But
    # we could refine the implementation in a few ways to get more meaningful
    # results (e.g. directly compute mmw for non-mmw cases, account for
    # metallicity in mmw calculations)

    # sanitize the arguments
    _out_dims = []
    ngrid_dims = 3 # number of dims used during construction

    def _sanitize_arr(arg, arg_name, grid_dim):
        ndim, size = np.ndim(arg), np.size(arg)
        if (ndim > 1) or (size == 0):
            raise ValueError(f"{arg_name!r} must be a scalar or a 1D array")
        elif ndim == 1:
            _out_dims.append(size)
        sanitized = np.atleast_1d(arg).astype(float)
        sanitized.shape = tuple(1 if i != grid_dim else -1
                                for i in range(ngrid_dims))
        return sanitized

    target_T_cgs = _sanitize_arr(target_T_cgs_arr, "target_T_cgs_arr", 0)
    target_Hmassdensity_cgs = _sanitize_arr(target_Hmassdensity_cgs_arr,
                                            "target_Hmassdensity_cgs_arr", 1)
    if (target_metallicity_arr is None):
        assert metal_cooling_flag in [0,1]
        # converting metal_cooling_flag to float happens to give proper default
        target_metallicity = _sanitize_arr(float(metal_cooling_flag), "-", 2)
    else:
        target_metallicity = _sanitize_arr(target_metallicity_arr,
                                           "target_metallicity_arr", 2)
    target_T_cgs, target_Hmassdensity_cgs, target_metalicity \
        = np.broadcast_arrays(target_T_cgs, target_Hmassdensity_cgs,
                              target_metallicity)

    tiny_number = 1e-20

    out = {}
    out["density"] = (target_Hmassdensity_cgs /
                      (nominal_HMassFrac * code_units['density_units']))

    if primordial_chem_flag > 0:
        out["HI"] = nominal_HMassFrac * out["density"]
        out["HII"] = tiny_number * out["density"]
        out["HeI"] = (1.0 - nominal_HMassFrac) * out["density"]
        out["HeII"] = tiny_number * out["density"]
        out["HeIII"] = tiny_number * out["density"]
        out["de"] = tiny_number * out["density"]
    if primordial_chem_flag > 1:
        out["H2I"] = tiny_number * out["density"]
        out["H2II"] = tiny_number * out["density"]
        out["HM"] = tiny_number * out["density"]
    if primordial_chem_flag > 2:
        out["DI"] = nominal_DeuteriumToHydrogenMassRatio * out["HI"]
        out["DII"] = tiny_number * out["density"]
        out["HDI"] = tiny_number * out["density"]
    if metal_cooling_flag == 1:
        out["metal"] = (out["density"] * nominal_SolarMetalFractionByMass *
                        target_metallicity)

    out["x-velocity"] = 0.0 * out["density"]
    out["y-velocity"] = 0.0 * out["density"]
    out["z-velocity"] = 0.0 * out["density"]

    initial_density_cgs = out["density"] * code_units['density_units']
    ndens_cgs_estimate = initial_density_cgs / (mass_hydrogen_cgs *
                                                nominal_initial_mmw)

    pressure_cgs = (
        ndens_cgs_estimate * boltzmann_constant_cgs * target_T_cgs)

    energy_cgs = (
        pressure_cgs / ((nominal_initial_gamma - 1.0) * initial_density_cgs))
    energy_units = (code_units["length_units"] / code_units["time_units"]) ** 2
    out["energy"] = energy_cgs / energy_units

    for k in out:
        out[k].shape = tuple(_out_dims)
    return out

