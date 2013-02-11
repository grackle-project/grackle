### Comment out Lines 535 to 544 of solve_rate_cool.F

import copy

def check_convergence(fc1, fc2, fields=None, tol=0.1):
    if fields is None:
        fields = ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                  "H2I", "H2II", "DI", "DII", "HDI", "de"]
    for field in fields:
        if (np.abs(fc1[field] - fc2[field]) / fc1[field] > tol).any():
            return False
    return True

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
my_chemistry.length_units = 1.0e16
my_chemistry.time_units = 3.15569e7 * 1e6
my_chemistry.a_units = 1.0

energy_units = (my_chemistry.length_units /
                my_chemistry.time_units)**2.0

a_value = 1.0

my_chemistry.initialize(a_value)

fc = FluidContainer(my_chemistry, 81)
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
fc["metal"][:] = 1.e-2 * fc["density"]

temperature_units = mass_h * ((my_chemistry.length_units/
                              my_chemistry.time_units)**2) / kboltz

initial_energy = np.logspace(1, 9, 81) / temperature_units

fc["energy"] = np.copy(initial_energy)
fc["x-velocity"][:] = 0.0
fc["y-velocity"][:] = 0.0
fc["z-velocity"][:] = 0.0

fc_last = copy.deepcopy(fc)

dt = 0.1
my_time = 0.0
while True:
    print "t = %e." % my_time
    for field in ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                  "H2I", "H2II", "DI", "DII", "HDI", "de"]:
        fc_last[field] = np.copy(fc[field])
    solve_chemistry(fc, a_value, dt)
    fc["energy"] = np.copy(initial_energy)    
    if check_convergence(fc, fc_last):
        break
    my_time += dt

calculate_temperature(fc)
calculate_cooling_time(fc, a_value, dt)
cooling_rate = fc["energy"] / fc["cooling_time"]

from matplotlib import pyplot
pyplot.loglog(fc['temperature'], cooling_rate)
pyplot.savefig('cooling_rate.png')
