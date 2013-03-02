import copy

def check_convergence(fc1, fc2, fields=None, tol=0.10):
    if fields is None:
        fields = ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                  "H2I", "H2II", "DI", "DII", "HDI", "de"]
    for field in fields:
        convergence = np.abs(fc1[field] - fc2[field]) / fc1[field]
        if (convergence > tol).any():
            print "Max change %s: %e." % (field, convergence.max())
            return False
    return True

kboltz      = 1.3806504e-16
mass_h      = 1.67262171e-24   
mass_e      = 9.10938215e-28
pi_val      = 3.14159265
hplanck     = 6.6260693e-27
ev2erg      = 1.60217653e-12
c_light     = 2.99792458e10
GravConst   = 6.6726e-8
sigma_sb    = 5.670373e-5
SolarMass   = 1.9891e33
Mpc         = 3.0857e24
kpc         = 3.0857e21
pc          = 3.0857e18
yr_to_s     = 3.15569e7

from pygrackle.grackle_wrapper import *
from pygrackle.fluid_container import FluidContainer

tiny_number = 1e-20

my_chemistry = chemistry_data()
my_chemistry.use_chemistry = 1
my_chemistry.with_radiative_cooling = 0
my_chemistry.primordial_chemistry = 3
my_chemistry.metal_cooling = 1
my_chemistry.grackle_data_file = "CloudyData_UVB=HM2012.h5"
my_chemistry.include_metal_heating = 1

my_chemistry.comoving_coordinates = 0
my_chemistry.density_units = mass_h
my_chemistry.length_units = 1.0e16
my_chemistry.time_units = yr_to_s * 1e6
my_chemistry.a_units = 1.0

energy_units = (my_chemistry.length_units /
                my_chemistry.time_units)**2.0

a_value = 1.0

my_chemistry.initialize(a_value)

n_points = 200
fc = FluidContainer(my_chemistry, n_points)
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
fc["metal"][:] = 0.02041 * fc["density"]

temperature_units = mass_h * ((my_chemistry.length_units/
                              my_chemistry.time_units)**2) / kboltz

initial_energy = np.logspace(1, 9, n_points) / temperature_units

fc["energy"] = np.copy(initial_energy)
fc["x-velocity"][:] = 0.0
fc["y-velocity"][:] = 0.0
fc["z-velocity"][:] = 0.0

fc_last = copy.deepcopy(fc)

dt = 0.1
my_time = 0.0
while True:
    print "t = %.2f Myr." % (my_time * my_chemistry.time_units /
                           (yr_to_s * 1e6))
    for field in ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                  "H2I", "H2II", "DI", "DII", "HDI", "de"]:
        fc_last[field] = np.copy(fc[field])
    solve_chemistry(fc, a_value, dt)
    fc["energy"] = np.copy(initial_energy)   
    if check_convergence(fc, fc_last):
        break
    my_time += dt

calculate_temperature(fc)
calculate_cooling_time(fc, a_value)
  
xbase1 = my_chemistry.length_units/(a_value*my_chemistry.a_units)
dbase1   = my_chemistry.density_units*(a_value*my_chemistry.a_units)**3
cool_unit = (my_chemistry.a_units**5 * xbase1**2 * mass_h**2) / \
  (my_chemistry.time_units**3 * dbase1)

cooling_rate = cool_unit * fc["density"] * fc["energy"] / \
  fc["cooling_time"]
  
from matplotlib import pyplot
pyplot.loglog(fc['temperature'], cooling_rate)
pyplot.xlabel('T [K]')
pyplot.ylabel('$\\Lambda$ [erg s$^{-1}$ cm$^{-3}$]')
pyplot.ylim(1e-30, 1e-21)
output_file = 'cooling_rate.png'
print "Writing %s." % output_file
pyplot.savefig(output_file)
