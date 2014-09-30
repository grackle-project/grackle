########################################################################
#
# Cooling rate example script
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this 
# software.
########################################################################

import copy
import numpy as np
import sys

from pygrackle.grackle_wrapper import *
from pygrackle.fluid_container import FluidContainer

from utilities.api import \
     setup_fluid_container, \
     set_cosmology_units, \
     get_cooling_units

from utilities.physical_constants import \
     mass_hydrogen_cgs, \
     sec_per_Gyr, \
     cm_per_mpc

# Set solver parameters
my_chemistry = chemistry_data()
my_chemistry.use_grackle = 1
my_chemistry.with_radiative_cooling = 0
my_chemistry.primordial_chemistry = 0
my_chemistry.metal_cooling = 1
my_chemistry.UVbackground = 1
my_chemistry.grackle_data_file = "CloudyData_UVB=HM2012.h5"
#my_chemistry.grackle_data_file = "CloudyData_noUVB.h5"

# Set units
my_chemistry.comoving_coordinates = 1 # proper units
my_chemistry.a_units = 1.0
my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
my_chemistry.time_units = sec_per_Gyr          # 1 Gyr in s
my_chemistry.velocity_units = my_chemistry.length_units / my_chemistry.time_units

current_redshift = 0.0

# Call convenience function for setting up a fluid container.
# This container holds the solver parameters, units, and fields.
temperature = np.logspace(1, 9, 200)
fc = setup_fluid_container(my_chemistry, 
                           temperature=temperature,
                           current_redshift=current_redshift,
                           converge=False)

a_value = 1.0 / (1.0 + current_redshift) / my_chemistry.a_units

calculate_temperature(fc)
calculate_cooling_time(fc, a_value)

cool_unit = get_cooling_units(my_chemistry, current_redshift)
density_proper = fc["density"] / \
  (my_chemistry.a_units * a_value)**(3*my_chemistry.comoving_coordinates)
cooling_rate = cool_unit * fc["energy"] / \
  np.abs(fc["cooling_time"]) / density_proper

t_sort = np.argsort(fc["temperature"])
  
from matplotlib import pyplot
pyplot.loglog(fc['temperature'][t_sort], cooling_rate[t_sort])
pyplot.xlabel('T [K]')
pyplot.ylabel('$\\Lambda$ [erg s$^{-1}$ cm$^{3}$]')
output_file = 'table_cooling_rate.png'
print "Writing %s." % output_file
pyplot.savefig(output_file)
