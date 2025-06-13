########################################################################
#
# Functions for evolving a fluid container using Grackle
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from collections import defaultdict
import numpy as np

from .physical_constants import \
    gravitational_constant_cgs, \
    sec_per_year

def evolve_freefall(fc, final_density, safety_factor=0.01,
                    include_pressure=True):
    my_chemistry = fc.chemistry_data

    # Set units of gravitational constant
    gravitational_constant = (
        4.0 * np.pi * gravitational_constant_cgs *
        my_chemistry.density_units * my_chemistry.time_units**2)

    # some constants for the analytical free-fall solution
    freefall_time_constant = np.power(((32. * gravitational_constant) /
                                       (3. * np.pi)), 0.5)

    data = defaultdict(list)
    current_time = 0.0
    while fc["density"][0] * my_chemistry.density_units < final_density:
        # calculate timestep based on free-fall solution
        dt = safety_factor * \
          np.power(((3. * np.pi) /
                    (32. * gravitational_constant *
                     fc["density"][0])), 0.5)

        add_to_data(fc, data, extra={"time": current_time})

        # compute the new density using the modified
        # free-fall collapse as per Omukai et al. (2005)
        if include_pressure:
            force_factor = calculate_collapse_factor(data["pressure"], data["density"])
        else:
            force_factor = 0.
        data["force_factor"].append(force_factor)

        # calculate new density from altered free-fall solution
        new_density = np.power((np.power(fc["density"][0], -0.5) -
                                (0.5 * freefall_time_constant * dt *
                                 np.power((1 - force_factor), 0.5))), -2.)

        print("Evolve Freefall - t: %e yr, rho: %e g/cm^3, T: %e K." %
              ((current_time * my_chemistry.time_units / sec_per_year),
               (fc["density"][0] * my_chemistry.density_units),
               fc["temperature"][0]))

        # use this to multiply by elemental densities if you are tracking those
        density_ratio = new_density / fc["density"][0]

        # update densities
        for field in fc.density_fields:
            fc[field] *= density_ratio

        # now update energy for adiabatic heating from collapse
        fc["internal_energy"][0] += (my_chemistry.Gamma - 1.) * \
          fc["internal_energy"][0] * freefall_time_constant * \
          np.power(fc["density"][0], 0.5) * dt

        fc.solve_chemistry(dt)

        # update time
        current_time += dt

    for field in data:
        data[field] = np.squeeze(np.array(data[field]))
    return fc.finalize_data(data=data)

def calculate_collapse_factor(pressure, density):
    # Calculate the effective adiabatic index, dlog(p)/dlog(rho).

    if len(pressure) < 3:
        return np.array([0.])

    # compute dlog(p) / dlog(rho) using last two timesteps
    gamma_eff = np.log10(pressure[-1] / pressure[-2]) / \
        np.log10(density[-1] / density[-2])

    # compute a higher order derivative if more than two points available
    if len(pressure) > 2:
        gamma_eff += 0.5 * ((np.log10(pressure[-2] / pressure[-3]) /
                             np.log10(density[-2] / density[-3])) - gamma_eff)

    gamma_eff = np.clip(gamma_eff, a_min=0, a_max=4/3)

    # Equation 9 of Omukai et al. (2005)
    if gamma_eff < 0.83:
        force_factor = gamma_eff * 0
    elif gamma_eff < 1.0:
        force_factor = 0.6 + 2.5 * (gamma_eff - 1) - \
            6.0 * np.power((gamma_eff - 1.0), 2.)
    else:
        force_factor = 1.0 + 0.2 * (gamma_eff - (4./3.)) - \
            2.9 * np.power((gamma_eff - (4./3.)), 2.)
    force_factor = np.clip(force_factor, a_min=0, a_max=0.95)
    return force_factor

def add_to_data(fc, data, extra=None):
    """
    Add current fluid container values to the data structure.
    """

    for field in fc.all_fields:
        if field not in fc.input_fields:
            func = getattr(fc, f"calculate_{field}")
            if func is None:
                raise RuntimeError(f"No function for calculating {field}.")
            func()
        data[field].append(fc[field].copy())

    if extra is not None:
        for field in extra:
            data[field].append(extra[field])
