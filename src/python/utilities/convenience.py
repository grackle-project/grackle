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
    for field in fields:
        if not field in fc1 or not field in fc2: continue
        convergence = np.abs(fc1[field] - fc2[field]) / fc1[field]
        if (convergence > tol).any():
            print "Max change %s: %.10e." % (field, convergence.max())
            return False
    return True

def setup_fluid_container(my_chemistry, 
                          density=mass_hydrogen_cgs,
                          current_redshift=0,
                          n_points=200, 
                          hydrogen_mass_fraction=0.76,
                          metal_mass_fraction=0.02041,
                          d_to_h_ratio=3.4e-5,
                          converge=False,
                          max_iterations=1000):
    """
    Initialize a fluid container with a constant density and smoothly 
    increasing temperature from 10 K to 1e9 K.  Optionally, iterate the 
    chemistry solver until the species fractions converge.  Return 
    the fluid container.
    """

    a_value = 1.0 / (1.0 + current_redshift) / my_chemistry.a_units

    my_value = my_chemistry.initialize(a_value)
    if not my_value:
        print "Error initializing chemistry."
        sys.exit(0)

    my_chemistry.update_UVbackground(a_value)

    tiny_number = 1e-20
    fc = FluidContainer(my_chemistry, n_points)
    fc["density"][:] = density / my_chemistry.density_units
    fc["HI"][:] = hydrogen_mass_fraction * fc["density"]
    fc["HII"][:] = tiny_number * fc["density"]
    fc["HeI"][:] = (1.0 - hydrogen_mass_fraction) * fc["density"]
    fc["HeII"][:] = tiny_number * fc["density"]
    fc["HeIII"][:] = tiny_number * fc["density"]
    fc["de"][:] = tiny_number * fc["density"]
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
    initial_energy = np.logspace(4, 9, n_points) / temperature_units
    fc["energy"] = np.copy(initial_energy)
    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    fc_last = fc.copy()

    dt = 0.01 * sec_per_Myr / my_chemistry.time_units
    my_time = 0.0
    i = 0
    while converge and i < max_iterations:
        print "t = %.2f Myr." % (my_time * my_chemistry.time_units /
                                 sec_per_Myr)
        for field in ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                      "H2I", "H2II", "DI", "DII", "HDI", "de"]:
            if field in fc:
                fc_last[field] = np.copy(fc[field])
        solve_chemistry(fc, a_value, dt)
        fc["energy"] = np.copy(initial_energy)   
        if check_convergence(fc, fc_last):
            break
        my_time += dt
        i += 1

    if i > max_iterations:
        print "ERROR: solver did not converge in %d iterations." % max_iterations
        
    return fc
