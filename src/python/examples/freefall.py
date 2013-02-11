# This is a translation of freefall.C from the clib distribution

kboltz      = 1.3806504e-16
mass_h      = 1.67262171e-24   
mass_e      = 9.10938215e-28
pi_val      = 3.14159265
hplanck     = 6.6260693e-27
ev2erg      = 1.60217653e-12
c_light     = 2.99792458e10
GravConst   = 6.67428e-8
sigma_sb    = 5.670373e-5
SolarMass   = 1.9891e33
Mpc         = 3.0857e24
kpc         = 3.0857e21
pc          = 3.0857e18

from pygrackle.grackle_wrapper import *
from pygrackle.fluid_container import FluidContainer

tiny_number = 1e-20

my_chemistry = chemistry_data()
my_chemistry.use_chemistry = 1
my_chemistry.primordial_chemistry = 3
my_chemistry.metal_cooling = 1
my_chemistry.cloudy_table_file = "solar_2008_3D_metals.h5";

my_chemistry.comoving_coordinates = 0
my_chemistry.density_units = 1.67e-24
my_chemistry.length_units = 1.0
my_chemistry.time_units = 1.0e12
my_chemistry.a_units = 1.0

energy_units = (my_chemistry.length_units /
                my_chemistry.time_units)**2.0

gravitational_constant = (4.0 * 3.1415926 * 6.6726e-8 * 
  my_chemistry.density_units * my_chemistry.time_units**2)

a_value = 1.0

my_chemistry.initialize(a_value)


my_chemistry.UVbackground = 1;
my_chemistry.UVbackground_file = "UVB_rates_HM2012.hdf5";
my_chemistry.initialize_UVbackground()
my_chemistry.update_UVbackground(a_value)


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


temperature_units = mass_h * ((my_chemistry.length_units/
                              my_chemistry.time_units)**2) / kboltz

fc["energy"][:] = 1000. / temperature_units
fc["x-velocity"][:] = 0.0
fc["y-velocity"][:] = 0.0
fc["z-velocity"][:] = 0.0

timestep_fraction = 0.1
current_time = 0.0
while fc["density"][0] < 1.e10:
    dt = timestep_fraction * \
      np.power(((3. * np.pi) / 
                (32. * gravitational_constant * fc["density"][0])), 0.5)

    calculate_temperature(fc)
    print fc["density"][0], fc["temperature"][0]
    
    density_ratio = np.power((freefall_constant - 
                              (0.5 * freefall_time_constant * current_time)), -2.) / \
                              fc["density"][0]
                              
    for field in fc.density_fields:
        fc[field] *= density_ratio

    fc["energy"][0] += (my_chemistry.Gamma - 1.) * fc["energy"][0] * \
      freefall_time_constant * np.power(fc["density"][0], 0.5) * dt

    solve_chemistry(fc, a_value, dt)

    current_time += dt
