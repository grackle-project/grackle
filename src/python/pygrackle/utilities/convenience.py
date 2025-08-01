########################################################################
#
# Python wrapper convenience functions
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np
import sys

from pygrackle.fluid_container import \
    FluidContainer

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr

def check_convergence(fc1, fc2, fields=None, tol=0.01):
    "Check for fields to be different by less than tol."

    if fields is None:
        fields = fc1.density_fields
    max_field = None
    max_val = 0.0
    for field in fields:
        if field not in fc2:
            continue
        convergence = np.max(np.abs(fc1[field] - fc2[field]) / fc1[field])
        if convergence > max_val:
            max_val = convergence
            max_field = field
    if np.any(max_val > tol):
        sys.stderr.write("max change - %5s: %.10e." % (max_field, max_val))
        return False
    return True

def setup_fluid_container(my_chemistry,
                          density=mass_hydrogen_cgs,
                          temperature=None,
                          state="ionized",
                          metal_mass_fraction=None,
                          dust_to_gas_ratio=None,
                          converge=False,
                          tolerance=0.01,
                          max_iterations=10000):
    """
    Initialize a fluid container using settings from a chemistry_data object.

    By default, initialize with a constant density and smoothly increasing
    temperature from 1e4 K to 1e9 K. Optionally, iterate the chemistry solver
    until the species fractions converge.

    The state of the gas can be set to "neutral" or "ionized" initialize the
    gas as fully neutral or fully ionized. Molecular fractions are always
    initialized to effectively zero. Note, if the primordial_chemistry
    parameter is set to 0, this option is ignored as only the total density
    and metal density are followed in this mode.

    Parameters
    ----------
    my_chemistry : chemistry_data
        Struct of Grackle runtime parameters.
    density : optional, float
        The density in CGS for all elements in the fluid container.
        Default: 1 hydrogen mass per cm^3 (i.e., ~1.67e-24 g/cm^3)
    temperature : optional, float or array of floats
        Temperature values in K.
    state : optional, string
        The initial ionization state of the gas. Either "neutral" to set
        all ionized species to effectively zero or "ionized" to set all
        neutral species to effectively zero.
    metal_mass_fraction : optional, float
        The mass fraction of gas in gas-phase metals.
        Default: 1e-20.
    dust_to_gas_ratio : optional, float
        The ratio of dust mass density to total gas density.
        Default: 1e-20.
    converge : optional, bool
        If True, iterate the solver until the chemical species reach
        equilibrium for the given temperatures.
        Default: False.
    tolerance : optional, float
        The maximum fractional change in a species density allowed to
        be considered converged.
        Default: 0.01.
    max_iterations : optional, int
        The maximum iterations to try for reaching convergence.
        Default: 10000.

    Returns
    -------
    fc : FluidContainer
        A fully initialized FluidContainer object.
    """

    rval = my_chemistry.initialize()
    if rval == 0:
        raise RuntimeError("Failed to initialize chemistry_data.")

    tiny_number = 1e-20
    if metal_mass_fraction is None:
        metal_mass_fraction = tiny_number
    if dust_to_gas_ratio is None:
        dust_to_gas_ratio = tiny_number
    if temperature is None:
        n_points = 200
        temperature = np.logspace(4, 9, n_points)
    else:
        if not isinstance(temperature, np.ndarray):
            temperature = np.array([temperature])
        n_points = temperature.size

    fc = FluidContainer(my_chemistry, n_points)
    fh = my_chemistry.HydrogenFractionByMass
    d2h = my_chemistry.DeuteriumToHydrogenRatio

    metal_free = 1 - metal_mass_fraction
    H_total = fh * metal_free
    He_total = (1 - fh) * metal_free
    # someday, maybe we'll include D in the total
    D_total = H_total * d2h

    fc_density = density / my_chemistry.density_units
    tiny_density = tiny_number * fc_density

    state_vals = {
        "density": fc_density,
        "metal_density": metal_mass_fraction * fc_density,
        "dust_density": dust_to_gas_ratio * fc_density
    }

    if state == "neutral":
        state_vals["HI_density"] = H_total * fc_density
        state_vals["HeI_density"] = He_total * fc_density
        state_vals["DI_density"] = D_total * fc_density
    elif state == "ionized":
        state_vals["HII_density"] = H_total * fc_density
        state_vals["HeIII_density"] = He_total * fc_density
        state_vals["DII_density"] = D_total * fc_density
        # ignore HeII since we'll set it to tiny
        state_vals["e_density"] = state_vals["HII_density"] + \
          state_vals["HeIII_density"] / 2
    else:
        raise ValueError("State must be either neutral or ionized.")

    for field in fc.density_fields:
        fc[field][:] = state_vals.get(field, tiny_density)

    fc.calculate_mean_molecular_weight()
    fc["internal_energy"] = temperature / \
        fc.chemistry_data.temperature_units / \
        fc["mean_molecular_weight"] / (my_chemistry.Gamma - 1.0)
    fc["x_velocity"][:] = 0.0
    fc["y_velocity"][:] = 0.0
    fc["z_velocity"][:] = 0.0

    fc_last = fc.copy()
    # disable cooling to iterate to equilibrium
    val = fc.chemistry_data.with_radiative_cooling
    fc.chemistry_data.with_radiative_cooling = 0

    my_time = 0.0
    i = 0
    while converge and i < max_iterations:
        fc.calculate_cooling_time()
        dt = 0.1 * np.abs(fc["cooling_time"]).min()
        sys.stderr.write("t: %.3f Myr, dt: %.3e Myr, " % \
                         ((my_time * my_chemistry.time_units / sec_per_Myr),
                          (dt * my_chemistry.time_units / sec_per_Myr)))
        for field in fc.density_fields:
            fc_last[field] = np.copy(fc[field])
        fc.solve_chemistry(dt)
        fc.calculate_mean_molecular_weight()
        fc["internal_energy"] = temperature / \
            fc.chemistry_data.temperature_units / fc["mean_molecular_weight"] / \
            (my_chemistry.Gamma - 1.0)
        converged = check_convergence(fc, fc_last, tol=tolerance)
        if converged:
            sys.stderr.write("\n")
            break
        sys.stderr.write("\r")
        my_time += dt
        i += 1

    fc.chemistry_data.with_radiative_cooling = val
    if i >= max_iterations:
        raise RuntimeError(
            f"ERROR: solver did not converge in {max_iterations} iterations.")

    return fc


def infer_eint_for_temperature_tabulated(
    chem,
    target_temperature,
    density,
    metal_density,
    *,
    cgs_density=False,
    xtol=None,
    rtol=None,
    maxiter=None
):
    """Infers internal energy assocaited with the target temperature

    Parameters
    ----------
    chem
        Chemistry object
    target_temperature
        Holds the target temperature (units of Kelvin)
    density
        Holds the total mass density (units are determined by the
        ``cgs_density`` kwarg)
    metal_density
        Holds the metal mass density (units are determined by the
        ``cgs_density`` kwarg)

    Note
    ----
    This function uses numerical root finding.

    While we could theoretically lookup a temperature by directly reading
    the cloudy table, that result wouldn't be exactly right. As the Grackle
    paper explains, Grackle uses a "dampener" when it computes the
    temperature a cloudy table. Consequently, it
    is "more correct" to use numerical root finding.
    """
    # TODO: don't require scipy (what we are doing is **REALLY** simple
    from scipy.optimize import root_scalar

    densityU = chem.density_units
    eintU = chem.get_velocity_units()**2
    gm1 = chem.Gamma - 1.0

    density = np.array(density)
    metal_density = np.array(metal_density)
    target_temperature = np.array(target_temperature)

    # todo: get a little more flexible about broadcasting
    # todo: allow metal_density to be omitted
    assert (density.ndim == 1) and (density.size > 0)
    assert density.shape == metal_density.shape == target_temperature.shape

    if chem.primordial_chemistry != 0:
        raise ValueError("Grackle wasn't configured in tabulated mode")

    # ensure that density & metal_density reference values in code units
    if cgs_density:
        density, metal_density = density / densityU, metal_density / densityU
    else:
        density, metal_density = density, metal_density

    # allocate the output array
    eint = np.empty(shape=density.shape, dtype="f8")

    # we are going to repeatedly call scipy.optimize.root_scalar. Let's set some stuff
    # up that we will reuse every time...
    fc = FluidContainer(chem, n_vals=1)

    def f(specific_internal_energy, target_T):
        fc["internal_energy"][0] = specific_internal_energy
        fc.calculate_temperature()
        return fc["temperature"] - target_T

    common_kwargs = {
        "f": f, "method": "bisect", "xtol": xtol, "rtol": rtol, "maxiter": maxiter
    }

    # move onto the actual loop:
    for i in range(density.size):
        # lets get bounds on the specific internal energy:
        # -> we could give more precise bounds on the mmw, and the way
        #    that metal density influences the value by reading the code
        #    or reading the discussion in GH-issue#67 (the paper actually
        #    discusses this, but actually had a mistake)
        mmw_bounds = np.array([0.6, 2.0])
        eint_bounds_cgs = (
            (boltzmann_constant_cgs * target_temperature[i]) /
            (gm1 * mmw_bounds * mass_hydrogen_cgs)
        )
        # convert to code units
        eint_bounds = eint_bounds_cgs / eintU
 
        # define the function used in root finding
        fc["density"][0] = density[i]
        fc["metal_density"][0] = metal_density[i]
 
        root_result = root_scalar(
            args=(target_temperature[i],), bracket=eint_bounds, **common_kwargs
        )
        if root_result.converged:
            eint[i] = root_result.root
        else:
            eint[i] = np.nan
    return eint
