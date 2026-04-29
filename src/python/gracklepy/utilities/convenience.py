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
import re
import roman
import sys

from gracklepy.fluid_container import \
    _element_masses, \
    FluidContainer

from gracklepy.utilities.atomic import \
    atomic_number, \
    mass_number, \
    primordial_elements, \
    solar_abundance
from gracklepy.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr

solar_mass_abundance = {element: solar_abundance[element] *
                        mass_number[element]
                        for element in solar_abundance
                        if element not in primordial_elements}
solar_metal_mass = sum(solar_mass_abundance.values())

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

def _get_appropriate_ion_field(fc, element, state):
    """
    Return the most ionized atomic species available for an element.
    """

    # this is entirely for deuterium, which we are ignoring
    if element not in atomic_number:
        return None

    if state == "neutral":
        order = range(atomic_number[element] + 2)
    elif state == "ionized":
        order = range(atomic_number[element] + 1, 0, -1)
    else:
        raise ValueError("State must be either neutral or ionized.")

    for i in order:
        fname = f"{element}{roman.toRoman(i)}_density"
        if fname in fc.density_fields:
            return fname

    # If we are tracking a field to represent the entire element, return that.
    fname = f"{element}_density"
    if fname in fc.density_fields:
        return fname

    raise ValueError(f"Element {element} not in fluid container")

def _setup_dust_densities(fc, state_vals, dust_to_gas_ratio):
    """
    Initialize the dust densities based on either the dust-to-gas ratio
    or on the injection metals.
    """

    if fc.chemistry_data.metal_chemistry == 0 or fc.chemistry_data.dust_species == 0:
        state_vals["dust_density"] = dust_to_gas_ratio * state_vals["density"]

    elif fc.chemistry_data.metal_chemistry == 1:
        metal_field = fc.inject_pathway_density_yield_fields[0]
        metal_density = state_vals[metal_field]
        dust_density = 0
        for gr, fmass in fc.chemistry_data._experimental_grain_inj_path_yields().items():
            state_vals[gr] = fmass[0] * metal_density
            dust_density += state_vals[gr]

        state_vals["dust_density"] = dust_density

    else:
        raise ValueError("metal_chemistry must be either 0 or 1.")


def _setup_metal_nuclide_densities(fc, state_vals, nuclide_densities):
    """
    Initialize the abundances of the metals based on the abundance pattern
    of the metal field. We either use a solar abundance pattern or that of
    an injection metal field if it is present and specifically selected.

    Note, here we are not setting state_vals. Instead, we are populating
    the nuclide_densities dict with will be used to set state_vals in
    setup_ion_fields.
    """

    # assume a solar abundance pattern
    if fc.chemistry_data.metal_chemistry == 0:
        metal_nuclide_fractions = {el: solar_mass_abundance[el] / solar_metal_mass
                                   for el in fc.elements if el not in primordial_elements}

    # use the injection metal data
    elif fc.chemistry_data.metal_chemistry == 1:
        metal_field = fc.inject_pathway_density_yield_fields[0]
        metal_nuclide_fractions = {el: fmass[0] for el, fmass in
          fc.chemistry_data._experimental_nuclide_gas_inj_path_yields().items()}

    else:
        raise ValueError("metal_chemistry must be either 0 or 1.")

    for el, fmass in metal_nuclide_fractions.items():
        if el not in fc.elements:
            continue
        nuclide_densities[el] = state_vals["metal_density"] * fmass

def _setup_ion_fields(fc, state_vals, nuclide_densities, state):
    """
    Initialize density fields for the ions that will the dominant
    species for either a neutral or ionized state.
    """

    state_vals["e_density"] = 0
    for el in nuclide_densities:
        fname = _get_appropriate_ion_field(fc, el, state)
        state_vals[fname] = nuclide_densities[el]

        # add to electron density
        if state == "ionized":
            reg = re.search(f"{el}(\w+)_density", fname)
            # skip density fields representing the whole element
            if reg is None:
                continue

            ion = reg.groups()[0]
            charge = roman.fromRoman(ion) - 1
            state_vals["e_density"] += nuclide_densities[el] * charge / \
              mass_number[el]

def _setup_inj_pathway_fields(state_vals: dict[str, float],
                              inj_pathway_yield_field_names: list[str]):
    """
    Helper function that updates state_vals to hold entries corresponding
    to the field(s) that specify the metals yielded by each injection pathway
    """
    n_pathways = len(inj_pathway_yield_field_names)
    if n_pathways == 0:
        return # nothing to be done
    elif n_pathways == 1:
        # put all the metal into the single yield we are following
        primary_pathway_yield_field = inj_pathway_yield_field_names[0]
    else:
        # we are following all possible metal yields, but for now
        # just put everything in the local ISM field.
        primary_pathway_yield_field = "local_ISM_metal_density"
        if primary_pathway_yield_field not in inj_pathway_yield_field_names:
            # at the time of writing, this shouldn't happen. But it will become
            # possible in the near future.
            raise RuntimeError(
                "We don't yet support the case where we have multiple injection "
                f"pathways but not the '{primary_pathway_yield_field}' field"
            )
    state_vals[primary_pathway_yield_field] = state_vals["metal_density"]


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

    fc_density = density / my_chemistry.density_units
    tiny_density = tiny_number * fc_density
    metal_free = 1 - metal_mass_fraction

    state_vals = {
        "density": fc_density,
        "metal_density": metal_mass_fraction * fc_density,
    }

    nuclide_densities = {
        "H": fh * metal_free,
        "He": (1 - fh) * metal_free,
    }
    # someday, maybe we'll include D in the total
    nuclide_densities["D"] = nuclide_densities["H"] * d2h
    for el in nuclide_densities:
        nuclide_densities[el] *= state_vals["density"]

    _setup_inj_pathway_fields(state_vals, fc.inject_pathway_density_yield_fields)
    _setup_metal_nuclide_densities(fc, state_vals, nuclide_densities)
    _setup_dust_densities(fc, state_vals, dust_to_gas_ratio)
    _setup_ion_fields(fc, state_vals, nuclide_densities, state)

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
