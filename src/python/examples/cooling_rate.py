import copy
import numpy as np
import sys

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

def check_convergence(fc1, fc2, fields=None, tol=0.01):
    if fields is None:
        fields = ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                  "H2I", "H2II", "DI", "DII", "HDI", "de"]
    for field in fields:
        convergence = np.abs(fc1[field] - fc2[field]) / fc1[field]
        if (convergence > tol).any():
            print "Max change %s: %.10e." % (field, convergence.max())
            return False
    return True

def set_cosmology_units(my_chemistry, hubble_constant=0.704,
                        omega_matter=0.268, omega_lambda=0.732,
                        current_redshift=0.0, initial_redshift=0.0,
                        comoving_box_size=1.0):
    """
    Greg Bryan's note on cosmology units:
    
    time:        utim = 1 / sqrt(4 * \pi * G * \rho_0 * (1+zri)^3)
    density:     urho = \rho_0 * (1+z)^3
    length:      uxyz = (1 Mpc) * box / h / (1+z)
    velocity:    uvel = uaye * uxyz / utim  (since u = a * dx/dt)
    (*)  temperature: utem = m_H * \mu / k * uvel**2
    a(t):        uaye = 1 / (1 + zri)
    
    where:
    box     - size of simulation box in Mpc/h
    zri     - initial redshift (start of simulation)
    \rho_0  = 3*\Omega_0*H_0^2/(8*\pi*G)
    Omega_0 - the fraction of non-relativistic matter at z=0
    
    Note that two definitions are dependent on redshift (urho
    and uxyz) so make sure to call this routine immediately
    before writing.
    
    * - the utem given below assumes that \mu = 1, so you must
    multiply the resulting temperature field by \mu.
    """
    
    my_chemistry.a_units = 1. / (1. + initial_redshift)
    my_chemistry.density_units = 1.8788e-29 * omega_matter * \
      np.power(hubble_constant,2) * np.power(1 + current_redshift,3)
    my_chemistry.length_units = 3.085678e24 * comoving_box_size / \
      hubble_constant / (1. + current_redshift)
    my_chemistry.time_units = 2.519445e17 / np.sqrt(omega_matter) / \
     hubble_constant / np.power(1 + initial_redshift, 1.5)
     
    velocity_units = my_chemistry.a_units * my_chemistry.length_units / \
      my_chemistry.time_units
    temperature_units = mass_h * velocity_units**2 / kboltz

    return temperature_units

from pygrackle.grackle_wrapper import *
from pygrackle.fluid_container import FluidContainer

tiny_number = 1e-20

my_chemistry = chemistry_data()
my_chemistry.use_chemistry = 1
my_chemistry.with_radiative_cooling = 0
my_chemistry.primordial_chemistry = 3
my_chemistry.metal_cooling = 1
my_chemistry.UVbackground = 1
my_chemistry.include_metal_heating = 1
my_chemistry.grackle_data_file = "CloudyData_UVB=HM2012.h5"
#my_chemistry.grackle_data_file = "CloudyData_noUVB.h5"

my_chemistry.comoving_coordinates = 1
initial_redshift = 99.0
current_redshift = 0.0
temperature_units = set_cosmology_units(my_chemistry, 
                                        current_redshift=current_redshift,
                                        initial_redshift=initial_redshift)
a_value = 1.0 / (1.0 + current_redshift) / my_chemistry.a_units

my_value = my_chemistry.initialize(a_value)
if not my_value:
    print "Error initializing chemistry."
    sys.exit(0)

my_chemistry.update_UVbackground(a_value)

n_points = 200
fc = FluidContainer(my_chemistry, n_points)
fc["density"][:] = 1.0 * mass_h / my_chemistry.density_units
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

initial_energy = np.logspace(1, 9, n_points) / temperature_units
fc["energy"] = np.copy(initial_energy)
fc["x-velocity"][:] = 0.0
fc["y-velocity"][:] = 0.0
fc["z-velocity"][:] = 0.0

fc_last = copy.deepcopy(fc)

dt = 1e4 * yr_to_s / my_chemistry.time_units
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
  
cooling_rate = cool_unit * fc["energy"] / \
  fc["cooling_time"] / (fc["density"] / a_value**3)

t_sort = np.argsort(fc["temperature"])
  
from matplotlib import pyplot
pyplot.loglog(fc['temperature'][t_sort], cooling_rate[t_sort])
pyplot.xlabel('T [K]')
pyplot.ylabel('$\\Lambda$ [erg s$^{-1}$ cm$^{3}$]')
output_file = 'cooling_rate.png'
print "Writing %s." % output_file
pyplot.savefig(output_file)
