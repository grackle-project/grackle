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
        fields = ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                  "H2I", "H2II", "DI", "DII", "HDI", "de"]
    max_field = None
    max_val = 0.0
    for field in fields:
        if field not in fc1 or field not in fc2:
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
                          hydrogen_mass_fraction=0.76,
                          metal_mass_fraction=0.02041,
                          d_to_h_ratio=3.4e-5,
                          converge=False, tolerance=0.01,
                          max_iterations=10000):
    """
    Initialize a fluid container with a constant density and smoothly
    increasing temperature from 10 K to 1e9 K.  Optionally, iterate the
    chemistry solver until the species fractions converge.  Return
    the fluid container.
    """

    rval = my_chemistry.initialize()
    if rval == 0:
        raise RuntimeError("Failed to initialize chemistry_data.")

    tiny_number = 1e-20
    if temperature is None:
        n_points = 200
        temperature = np.logspace(4, 9, n_points)
    else:
        n_points = temperature.size
    fc = FluidContainer(my_chemistry, n_points)
    fc["density"][:] = density / my_chemistry.density_units
    if my_chemistry.primordial_chemistry > 0:
        fc["HII"][:] = hydrogen_mass_fraction * fc["density"]
        fc["HI"][:] = tiny_number * fc["density"]
        fc["HeI"][:] = (1.0 - hydrogen_mass_fraction) * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
        fc["de"][:] = fc["HII"] + fc["HeII"] / 4.0 + fc["HeIII"] / 2.0
    if my_chemistry.primordial_chemistry > 1:
        fc["HM"][:] = tiny_number * fc["density"]
        fc["H2I"][:] = tiny_number * fc["density"]
        fc["H2II"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * d_to_h_ratio * fc["density"]
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    fc["metal"][:] = metal_mass_fraction * fc["density"]

    fc["energy"] = temperature / \
        fc.chemistry_data.temperature_units / \
        fc.calculate_mean_molecular_weight() / \
        (my_chemistry.Gamma - 1.0)
    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    fc_last = fc.copy()

    my_time = 0.0
    i = 0
    while converge and i < max_iterations:
        fc.calculate_cooling_time()
        dt = 0.1 * np.abs(fc["cooling_time"]).min()
        sys.stderr.write("t: %.3f Myr, dt: %.3e Myr, " % \
                         ((my_time * my_chemistry.time_units / sec_per_Myr),
                          (dt * my_chemistry.time_units / sec_per_Myr)))
        for field in ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                      "H2I", "H2II", "DI", "DII", "HDI", "de"]:
            if field in fc:
                fc_last[field] = np.copy(fc[field])
        fc.solve_chemistry(dt)
        mu = fc.calculate_mean_molecular_weight()
        fc["energy"] = temperature / \
            fc.chemistry_data.temperature_units / mu / \
            (my_chemistry.Gamma - 1.0)
        converged = check_convergence(fc, fc_last, tol=tolerance)
        if converged:
            sys.stderr.write("\n")
            break
        sys.stderr.write("\r")
        my_time += dt
        i += 1

    if i >= max_iterations:
        sys.stderr.write("ERROR: solver did not converge in %d iterations.\n" %
                         max_iterations)
        return None

    return fc
