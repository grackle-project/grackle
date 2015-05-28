########################################################################
#
# Cooling cell example script
#
#  This will initialize a single cell at a given temperature,  
#  iterate the cooling solver for a fixed time, and output the 
#  temperature vs. time.
#
#
# Copyright (c) 2015, Enzo/Grackle Development Team.
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
     get_temperature_units

from utilities.physical_constants import \
     mass_hydrogen_cgs, \
     sec_per_Myr, \
     sec_per_year, \
     cm_per_mpc

tiny_number = 1e-20

# Set initial values
density             = 0.1 # g /cm^3
initial_temperature = 1.e6 # K
final_time          = 100.0 # Myr

# Set solver parameters
my_chemistry = chemistry_data()
my_chemistry.use_grackle = 1
my_chemistry.with_radiative_cooling = 1
my_chemistry.primordial_chemistry = 0
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
my_chemistry.velocity_units = my_chemistry.length_units / \
  my_chemistry.time_units

a_value = 1.0

my_value = my_chemistry.initialize(a_value)
if not my_value:
    print "Error initializing chemistry."
    sys.exit(0)

fc = FluidContainer(my_chemistry, 1)
fc["density"][:] = density
if my_chemistry.primordial_chemistry > 0:
    fc["HI"][:] = 0.76 * fc["density"]
    fc["HII"][:] = tiny_number * fc["density"]
    fc["HeI"][:] = (1.0 - 0.76) * fc["density"]
    fc["HeII"][:] = tiny_number * fc["density"]
    fc["HeIII"][:] = tiny_number * fc["density"]
if my_chemistry.primordial_chemistry > 1:
    fc["H2I"][:] = tiny_number * fc["density"]
    fc["H2II"][:] = tiny_number * fc["density"]
    fc["HM"][:] = tiny_number * fc["density"]
    fc["de"][:] = tiny_number * fc["density"]
if my_chemistry.primordial_chemistry > 2:
    fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
    fc["DII"][:] = tiny_number * fc["density"]
    fc["HDI"][:] = tiny_number * fc["density"]
if my_chemistry.metal_cooling == 1:
    fc["metal"][:] = 0.1 * fc["density"] * \
      my_chemistry.SolarMetalFractionByMass
    
fc["x-velocity"][:] = 0.0
fc["y-velocity"][:] = 0.0
fc["z-velocity"][:] = 0.0

temperature_units = get_temperature_units(my_chemistry)
fc["energy"][:] = initial_temperature / temperature_units
calculate_temperature(fc, a_value)
fc["energy"][:] *= initial_temperature / fc["temperature"]

temperature_values = []
mu_values = []
time_values = []

timestep_fraction = 0.01
current_time = 0.0

while current_time < final_time:
    calculate_cooling_time(fc, a_value)
    dt = timestep_fraction * np.abs(fc["cooling_time"][0])

    calculate_temperature(fc, a_value)

    temperature_values.append(fc["temperature"][0])
    time_values.append(current_time * my_chemistry.time_units)
    mu_values.append(fc["temperature"][0] /
                     (fc["energy"][0] * (my_chemistry.Gamma - 1.) *
                      temperature_units))

    print "t: %e yr, T: %e K, dt: %e yr" % \
      ((current_time * my_chemistry.time_units / sec_per_year),
       fc["temperature"][0],
       (dt * my_chemistry.time_units / sec_per_year))

    solve_chemistry(fc, a_value, dt)
            
    current_time += dt

from matplotlib import pyplot

p1, = pyplot.loglog(time_values, temperature_values, 
                    color="k", label="T")
pyplot.xlabel("Time [s]")
pyplot.ylabel("T [K]")
output_file = "cooling_cell.png"

pyplot.twinx()
p2, = pyplot.semilogx(time_values, mu_values,
                      color="r", label="$\\mu$")
pyplot.ylabel("$\\mu$")
pyplot.legend([p1,p2],["T","$\\mu$"], fancybox=True,
              loc="center left")

print "Writing %s." % output_file
pyplot.savefig(output_file)

fh = file("cooling.out", "w")
fh.write("# t [s]      T[K]\n")
for t, T in zip(time_values, temperature_values):
    fh.write("%e %e\n" % (t, T))
fh.close()
