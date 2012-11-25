# This is a translation of freefall.C from the clib distribution

from pygrackle.grackle_wrapper import *
from pygrackle.fluid_container import FluidContainer

my_chemistry = chemistry_data()
my_chemistry.use_chemistry = 1
my_chemistry.primordial_chemistry = 2
#my_chemistry.metal_cooling = 1
#my_chemistry.cloudy_table_file = "solar_2008_3D_metals.h5";

my_chemistry.comoving_coordinates = 0
my_chemistry.density_units = 1.67e-24
my_chemistry.length_units = 1.0
my_chemistry.time_units = 1.0e12
my_chemistry.a_units = 1.0

energy_units = (my_chemistry.length_units / my_chemistry.time_units)**2.0

gravitational_constant = (4.0 * 3.1415926 * 6.6726e-8 * 
  my_chemistry.density_units * my_chemistry.time_units**2)

a_value = 1.0/(1.0+10.0)

my_chemistry.initialize(a_value)

fc = FluidContainer(my_chemistry, 64)
