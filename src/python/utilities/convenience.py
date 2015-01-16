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

from pygrackle.grackle_wrapper import *
from pygrackle.fluid_container import FluidContainer

from .units import \
     set_cosmology_units, \
     get_temperature_units
     
from utilities.physical_constants import \
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
        if not field in fc1 or not field in fc2: continue
        convergence = np.max(np.abs(fc1[field] - fc2[field]) / fc1[field])
        if convergence > max_val:
            max_val = convergence
            max_field = field
    if (max_val > tol).any():
        print "Max change %s: %.10e." % (max_field, max_val)
        return False
    return True

def setup_fluid_container(my_chemistry,
                          density=mass_hydrogen_cgs,
                          temperature=None,
                          current_redshift=0,
                          hydrogen_mass_fraction=0.76,
                          metal_mass_fraction=0.02041,
                          d_to_h_ratio=3.4e-5,
                          converge=False, tolerance=0.01,
                          max_iterations=10000, dt=None):
    """
    Initialize a fluid container with a constant density and smoothly 
    increasing temperature from 10 K to 1e9 K.  Optionally, iterate the 
    chemistry solver until the species fractions converge.  Return 
    the fluid container.
    """

    a_value = 1.0 / (1.0 + current_redshift) / my_chemistry.a_units

    my_value = my_chemistry.initialize(a_value)
    if not my_value:
        sys.stderr.write("Error initializing chemistry.\n")
        return None

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

    temperature_units = get_temperature_units(my_chemistry)
    fc["energy"] = temperature / temperature_units / \
      calculate_mean_molecular_weight(my_chemistry, fc)
    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    fc_last = fc.copy()

    if dt is None:
        dt = 0.01 * sec_per_Myr / my_chemistry.time_units
    my_time = 0.0
    i = 0
    while converge and i < max_iterations:
        print "t = %.3f Myr, dt = %.3e Myr" % \
          ((my_time * my_chemistry.time_units / sec_per_Myr),
           (dt * my_chemistry.time_units / sec_per_Myr))
        for field in ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                      "H2I", "H2II", "DI", "DII", "HDI", "de"]:
            if field in fc:
                fc_last[field] = np.copy(fc[field])
        solve_chemistry(fc, a_value, dt)
        mu = calculate_mean_molecular_weight(my_chemistry, fc)
        fc["energy"] = temperature / temperature_units / mu
        converged = check_convergence(fc, fc_last, tol=tolerance)
        if converged: break
        my_time += dt
        i += 1

    if i >= max_iterations:
        sys.stderr.write("ERROR: solver did not converge in %d iterations.\n" % 
                         max_iterations)
        return None
        
    return fc

def calculate_mean_molecular_weight(my_chemistry, fc):
    mu_metal = 16.0
    if my_chemistry.primordial_chemistry == 0:
        return calc_mu_table(fc['temperature'])
    mu = fc["HI"] + fc["HII"] + fc["de"] + \
      (fc["HeI"] + fc["HeII"] + fc["HeIII"]) / 4.
    if my_chemistry.primordial_chemistry > 1:
        mu += fc["HM"] + fc["H2I"] + fc["H2II"]
    if my_chemistry.primordial_chemistry > 2:
        mu += (fc["DI"] + fc["DII"]) / 2. + fc["HDI"] / 3.

    if my_chemistry.metal_cooling == 1:
        mu += fc["metal"] / mu_metal

    return fc["density"] / mu

def calculate_hydrogen_number_density(my_chemistry, fc):
    if my_chemistry.primordial_chemistry == 0:
        return my_chemistry.HydrogenFractionByMass * \
          fc["density"] * my_chemistry.density_units / mass_hydrogen_cgs
    nH = fc["HI"] + fc["HII"]
    if my_chemistry.primordial_chemistry > 1:
        nH += fc["HM"] + fc["H2I"] + fc["H2II"]
    if my_chemistry.primordial_chemistry > 2:
        nH += fc["HDI"] / 2.
    return nH * my_chemistry.density_units / mass_hydrogen_cgs


# ------------------------------------------------------------------
#   Calculate a tabulated approximation to mean molecular weight.
# ------------------------------------------------------------------
def calc_mu_table(temperature):
    tt = np.array([1.0e+01, 1.0e+02, 1.0e+03, 1.0e+04, 1.3e+04, 2.1e+04, \
                   3.4e+04, 6.3e+04, 1.0e+05, 1.0e+09])
    mt = np.array([1.18701555, 1.15484424, \
                   1.09603514, 0.9981496, 0.96346395, 0.65175895, \
                   0.6142901,  0.6056833, 0.5897776,  0.58822635])
    
    logttt= np.log(temperature)
    logmu = np.interp(logttt,np.log(tt),np.log(mt)) # linear interpolation in log-log space
    return np.exp(logmu)
