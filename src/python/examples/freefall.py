########################################################################
#
# Free-fall example script
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

from pygrackle.grackle_wrapper import \
     calculate_cooling_time, \
     calculate_temperature, \
     chemistry_data, \
     solve_chemistry

from pygrackle.fluid_container import FluidContainer

from utilities.api import \
     get_cooling_units, \
     get_temperature_units

from utilities.physical_constants import \
     gravitational_constant_cgs, \
     mass_hydrogen_cgs, \
     sec_per_Myr, \
     sec_per_year, \
     cm_per_mpc

tiny_number = 1e-20
     
# Set solver parameters
my_chemistry = chemistry_data()
my_chemistry.use_grackle = 1
my_chemistry.with_radiative_cooling = 1
my_chemistry.primordial_chemistry = 3
my_chemistry.metal_cooling = 1
my_chemistry.UVbackground = 1
my_chemistry.grackle_data_file = "CloudyData_UVB=HM2012.h5"
#my_chemistry.grackle_data_file = "CloudyData_noUVB.h5"

# Set units
my_chemistry.comoving_coordinates = 0 # proper units
my_chemistry.a_units = 1.0
my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
my_chemistry.time_units = sec_per_Myr          # 1 Myr in s
my_chemistry.velocity_units = my_chemistry.length_units / my_chemistry.time_units

# Set units of gravitational constant
gravitational_constant = (4.0 * np.pi * gravitational_constant_cgs * 
  my_chemistry.density_units * my_chemistry.time_units**2)

a_value = 1.0

my_value = my_chemistry.initialize(a_value)
if not my_value:
    print "Error initializing chemistry."
    sys.exit(0)

fc = FluidContainer(my_chemistry, 1)
fc["density"][:] = 1.0
fc["HI"][:] = 0.76 * fc["density"]
fc["HII"][:] = tiny_number * fc["density"]
fc["HM"][:] = tiny_number * fc["density"]
fc["HeI"][:] = (1.0 - 0.76) * fc["density"]
fc["HeII"][:] = tiny_number * fc["density"]
fc["HeIII"][:] = tiny_number * fc["density"]
fc["H2I"][:] = tiny_number * fc["density"]
fc["H2II"][:] = tiny_number * fc["density"]
fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
fc["DII"][:] = tiny_number * fc["density"]
fc["HDI"][:] = tiny_number * fc["density"]
fc["de"][:] = tiny_number * fc["density"]
fc["metal"][:] = 1.e-5 * fc["density"]

freefall_constant = np.power(fc["density"][0], -0.5)
freefall_time_constant = np.power(((32. * gravitational_constant) / 
                                   (3. * np.pi)), 0.5)


temperature_units = get_temperature_units(my_chemistry)

fc["energy"][:] = 1000. / temperature_units
fc["x-velocity"][:] = 0.0
fc["y-velocity"][:] = 0.0
fc["z-velocity"][:] = 0.0

density_values = []
temperature_values = []
time_values = []

timestep_fraction = 0.1
current_time = 0.0
while fc["density"][0] < 1.e10:
    dt = timestep_fraction * \
      np.power(((3. * np.pi) / 
                (32. * gravitational_constant * 
                 fc["density"][0])), 0.5)

    calculate_temperature(fc, a_value)

    density_values.append(fc["density"][0] * my_chemistry.density_units)
    temperature_values.append(fc["temperature"][0])
    time_values.append(current_time * my_chemistry.time_units / sec_per_year)
    
    print "t: %e yr, rho: %e g/cm^3, T: %e [K]" % \
      ((current_time * my_chemistry.time_units / sec_per_year),
       (fc["density"][0] * my_chemistry.density_units),
       fc["temperature"][0])
    
    density_ratio = np.power((freefall_constant - 
                              (0.5 * freefall_time_constant * 
                               current_time)), -2.) / \
                              fc["density"][0]
                              
    for field in fc.density_fields:
        fc[field] *= density_ratio

    fc["energy"][0] += (my_chemistry.Gamma - 1.) * fc["energy"][0] * \
      freefall_time_constant * np.power(fc["density"][0], 0.5) * dt

    solve_chemistry(fc, a_value, dt)
    
    current_time += dt

from matplotlib import pyplot

p1, = pyplot.loglog(time_values, density_values)
pyplot.xlabel('Time [yr]')
pyplot.ylabel('$\\rho$ [g cm$^{-3}$]')
pyplot.axis([1e6,4e7,1e-24,1e-13])

pyplot.twinx()
p2, = pyplot.loglog(time_values, temperature_values, color='r')
pyplot.ylabel('T [K]')
pyplot.axis([1e6,4e7,50.0,5e3])

pyplot.legend([p1,p2],['density','temperature'],fancybox=True,loc='center left')

output_file = 'freefall.png'
print "Writing %s." % output_file
pyplot.savefig(output_file)
