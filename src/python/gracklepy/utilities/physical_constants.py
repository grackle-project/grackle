#
# Physical Constants and Units Conversion Factors
#
# Values for these constants, unless otherwise noted, are drawn from IAU,
# IUPAC, and NIST data, whichever is newer.
# http://maia.usno.navy.mil/NSFA/IAU2009_consts.html
# http://goldbook.iupac.org/list_goldbook_phys_constants_defs.html
# http://physics.nist.gov/cuu/Constants/index.html

# Masses
mass_electron_cgs = 9.109382e-28  # g
amu_cgs           = 1.660538921e-24  # g
mass_hydrogen_cgs = 1.007947*amu_cgs  # g
mass_sun_cgs = 1.98841586e33  # g
# Velocities
speed_of_light_cgs = 2.99792458e10  # cm/s, exact

# Cross Sections
# 8*pi/3 (alpha*hbar*c/(2*pi))**2
cross_section_thompson_cgs = 6.65245854533e-25  # cm^2

# Charge
charge_proton_cgs = 4.8032056e-10  # esu = 1.602176487e-19  Coulombs

# Physical Constants
boltzmann_constant_cgs = 1.3806488e-16  # erg K^-1
gravitational_constant_cgs  = 6.67384e-8  # cm^3 g^-1 s^-2
planck_constant_cgs   = 6.62606957e-27  # erg s
stefan_boltzmann_constant_cgs = 5.670373e-5 # erg cm^-2 s^-1 K^-4
# The following value was calcualted assuming H = 100 km/s/Mpc.
# To get the correct value for your cosmological parameters,
# you'll need to multiply through by h^2
# [where h = H / (100 km/s/Mpc)].  See the Overdensity field in
# yt.data_objects.universal_fields.
rho_crit_now = 1.8788e-29  # g/cm^3 (cosmological critical density)


# Misc. Approximations
mass_mean_atomic_cosmology = 1.22
mass_mean_atomic_galactic = 2.3

# Conversion Factors:  X au * mpc_per_au = Y mpc
# length
mpc_per_mpc   = 1e0
mpc_per_kpc   = 1e-3
mpc_per_pc    = 1e-6
mpc_per_au    = 4.84813682e-12
mpc_per_rsun  = 2.253962e-14
mpc_per_miles = 5.21552871e-20
mpc_per_km    = 3.24077929e-20
mpc_per_cm    = 3.24077929e-25
kpc_per_cm    = mpc_per_cm / mpc_per_kpc
km_per_pc     = 3.08567758e13
km_per_m      = 1e-3
km_per_cm     = 1e-5
pc_per_cm     = 3.24077929e-19

m_per_fpc     = 0.0324077929

kpc_per_mpc   = 1.0 / mpc_per_kpc
pc_per_mpc    = 1.0 / mpc_per_pc
au_per_mpc    = 1.0 / mpc_per_au
rsun_per_mpc  = 1.0 / mpc_per_rsun
miles_per_mpc = 1.0 / mpc_per_miles
km_per_mpc    = 1.0 / mpc_per_km
cm_per_mpc    = 1.0 / mpc_per_cm
cm_per_kpc    = 1.0 / kpc_per_cm
cm_per_km     = 1.0 / km_per_cm
pc_per_km     = 1.0 / km_per_pc
cm_per_pc     = 1.0 / pc_per_cm

# time
# "IAU Style Manual" by G.A. Wilkins, Comm. 5, in IAU Transactions XXB (1989)
sec_per_Gyr  = 31.5576e15
sec_per_Myr  = 31.5576e12
sec_per_kyr  = 31.5576e9
sec_per_year = 31.5576e6
sec_per_day  = 86400.0
sec_per_hr   = 3600.0
day_per_year = 365.25

# temperature / energy
erg_per_eV = 1.602176562e-12
erg_per_keV = erg_per_eV * 1.0e3
K_per_keV = erg_per_keV / boltzmann_constant_cgs
keV_per_K = 1.0 / K_per_keV

#Short cuts
G = gravitational_constant_cgs
me = mass_electron_cgs
mp = mass_hydrogen_cgs
qp = charge_proton_cgs
mh = mp
clight = speed_of_light_cgs
kboltz = boltzmann_constant_cgs
hcgs = planck_constant_cgs
sigma_thompson = cross_section_thompson_cgs
Na = 1 / amu_cgs
