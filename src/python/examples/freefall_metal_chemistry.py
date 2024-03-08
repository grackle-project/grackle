########################################################################
#
# Free-fall example script
#
#
# Copyright (c) 2013-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import os
import yt

from pygrackle import \
    chemistry_data, \
    FluidContainer, \
    evolve_constant_density, \
    evolve_freefall

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    mass_electron_cgs, \
    sec_per_Myr, \
    cm_per_mpc

tiny_number = 1e-60

solar_abundance = {
    'H' : 1.00e+00, 'He': 1.00e-01, 'Li': 2.04e-09,
    'Be': 2.63e-11, 'B' : 6.17e-10, 'C' : 2.45e-04,
    'N' : 8.51e-05, 'O' : 4.90e-04, 'F' : 3.02e-08,
    'Ne': 1.00e-04, 'Na': 2.14e-06, 'Mg': 3.47e-05,
    'Al': 2.95e-06, 'Si': 3.47e-05, 'P' : 3.20e-07,
    'S' : 1.84e-05, 'Cl': 1.91e-07, 'Ar': 2.51e-06,
    'K' : 1.32e-07, 'Ca': 2.29e-06, 'Sc': 1.48e-09,
    'Ti': 1.05e-07, 'V' : 1.00e-08, 'Cr': 4.68e-07,
    'Mn': 2.88e-07, 'Fe': 2.82e-05, 'Co': 8.32e-08,
    'Ni': 1.78e-06, 'Cu': 1.62e-08, 'Zn': 3.98e-08}

atomic_mass = {
    'H' : 1.00794,   'He': 4.002602,  'Li': 6.941,
    'Be': 9.012182,  'B' : 10.811,    'C' : 12.0107,
    'N' : 14.0067,   'O' : 15.9994,   'F' : 18.9984032,
    'Ne': 20.1797,   'Na': 22.989770, 'Mg': 24.3050,
    'Al': 26.981538, 'Si': 28.0855,   'P' : 30.973761,
    'S' : 32.065,    'Cl': 35.453,    'Ar': 39.948,
    'K' : 39.0983,   'Ca': 40.078,    'Sc': 44.955910,
    'Ti': 47.867,    'V' : 50.9415,   'Cr': 51.9961,
    'Mn': 54.938049, 'Fe': 55.845,    'Co': 58.933200,
    'Ni': 58.6934,   'Cu': 63.546,    'Zn': 65.409}

atomic_number = {
    'H' : 1,  'He': 2,  'Li': 3,
    'Be': 4,  'B' : 5,  'C' : 6,
    'N' : 7,  'O' : 8,  'F' : 9,
    'Ne': 10, 'Na': 11, 'Mg': 12,
    'Al': 13, 'Si': 14, 'P' : 15,
    'S' : 16, 'Cl': 17, 'Ar': 18,
    'K' : 19, 'Ca': 20, 'Sc': 21,
    'Ti': 22, 'V' : 23, 'Cr': 24,
    'Mn': 25, 'Fe': 26, 'Co': 27,
    'Ni': 28, 'Cu': 29, 'Zn': 30}

solar_total_mass = sum(solar_abundance[a] * atomic_mass[a]
                    for a in solar_abundance)
solar_metal_mass = sum(solar_abundance[a] * atomic_mass[a]
                    for a in solar_abundance if a not in ["H", "He"])
solar_mass_fraction = \
  {a: solar_abundance[a] * atomic_mass[a] / solar_total_mass
   for a in solar_abundance}
primordial_mass_fraction = {"H": 0.76, "He": 0.24}

def get_mass_fraction(el, metallicity):
    if el not in ("H", "He"):
        return metallicity * solar_mass_fraction.get(el, 0)

    # For H/He interpolate as fraction of primordial mass
    # so that we always add up to 1
    XY_1 = 1 - solar_metal_mass / solar_total_mass
    XY_Z = 1 - metallicity * solar_metal_mass / solar_total_mass
    fXY_0 = primordial_mass_fraction.get(el, 0)
    fXY_1 = solar_mass_fraction.get(el, 0) / XY_1
    fXY_Z = (fXY_1 - fXY_0) * metallicity + fXY_0
    return fXY_Z * XY_Z

if __name__=="__main__":
    current_redshift = 0.

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    my_chemistry.primordial_chemistry = 4
    my_chemistry.UVbackground = 0
    my_chemistry.self_shielding_method = 0
    my_chemistry.H2_self_shielding = 0
    my_chemistry.Gamma = 5. / 3.
    my_chemistry.CaseBRecombination = 0
    my_chemistry.cie_cooling = 1
    my_chemistry.h2_optical_depth_approximation = 1
    my_chemistry.interstellar_radiation_field = 0.

    my_chemistry.metal_cooling = 1
    my_chemistry.metal_chemistry = 1
    my_dir = os.path.dirname(os.path.abspath(__file__))
    grackle_data_file = os.path.join(
        my_dir, "..", "..", "..", "input", "cloudy_metals_2008_3D.h5")
    my_chemistry.grackle_data_file = grackle_data_file
    my_chemistry.h2_on_dust = 1
    my_chemistry.use_dust_density_field = 1
    metallicity = 1e-4

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units  = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units   = cm_per_mpc        # 1 Mpc in cm
    my_chemistry.time_units     = sec_per_Myr       # 1 Myr in s
    my_chemistry.set_velocity_units()

    # set initial density and temperature
    initial_temperature = 50000. # start the gas at this temperature
    # then begin collapse
    initial_density     = 1.0e-1 * mass_hydrogen_cgs # g / cm^3
    # stopping condition
    final_density       = 1.e12 * mass_hydrogen_cgs

    rval = my_chemistry.initialize()

    fc = FluidContainer(my_chemistry, 1)
    fc["density"][:] = initial_density / my_chemistry.density_units
    fc["HI"][:] = get_mass_fraction("H", metallicity) * fc["density"]
    fc["HII"][:] = tiny_number * get_mass_fraction("H", metallicity) * fc["density"]
    fc["HeI"][:] = get_mass_fraction("He", metallicity) * fc["density"]
    fc["HeII"][:] = tiny_number * get_mass_fraction("He", metallicity) * fc["density"]
    fc["HeIII"][:] = tiny_number * get_mass_fraction("He", metallicity) * fc["density"]
    fc["de"][:] = 2e-4 * mass_electron_cgs / mass_hydrogen_cgs * fc["density"]
    if my_chemistry.primordial_chemistry > 1:
        fc["H2I"][:] = tiny_number * fc["density"]
        fc["H2II"][:] = tiny_number * fc["density"]
        fc["HM"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 3:
        fc["DM"][:] = tiny_number * fc["density"]
        fc["HDII"][:] = tiny_number * fc["density"]
        fc["HeHII"][:] = tiny_number * fc["density"]
    if my_chemistry.metal_cooling == 1:
        fc["metal"][:] = metallicity * fc["density"] * \
            my_chemistry.SolarMetalFractionByMass
    if my_chemistry.use_dust_density_field:
        fc["dust"][:] = metallicity * fc["density"] * \
            my_chemistry.local_dust_to_gas_ratio
    if my_chemistry.metal_chemistry > 0:
        # this is not exactly correct
        fc["CI"][:] = get_mass_fraction("C", metallicity) * fc["density"]
        fc["CII"][:] = tiny_number * fc["density"]
        fc["CO"][:] = tiny_number * fc["density"]
        fc["CO2"][:] = tiny_number * fc["density"]
        fc["OI"][:] = get_mass_fraction("O", metallicity) * fc["density"]
        fc["OH"][:] = tiny_number * fc["density"]
        fc["H2O"][:] = tiny_number * fc["density"]
        fc["O2"][:] = tiny_number * fc["density"]
        fc["SiI"][:] = get_mass_fraction("Si", metallicity) * fc["density"]
        fc["SiOI"][:] = tiny_number * fc["density"]
        fc["SiO2I"][:] = tiny_number * fc["density"]
        fc["CH"][:] = tiny_number * fc["density"]
        fc["CH2"][:] = tiny_number * fc["density"]
        fc["COII"][:] = tiny_number * fc["density"]
        fc["OII"][:] = tiny_number * fc["density"]
        fc["OHII"][:] = tiny_number * fc["density"]
        fc["H2OII"][:] = tiny_number * fc["density"]
        fc["H3OII"][:] = tiny_number * fc["density"]
        fc["O2II"][:] = tiny_number * fc["density"]

    fc["energy"][:] = initial_temperature / \
        fc.chemistry_data.temperature_units
    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    # timestepping safety factor
    safety_factor = 0.01

    # let the gas cool at constant density from the starting temperature
    # down to a lower temperature to get the species fractions in a
    # reasonable state.
    cooling_temperature = 100.
    data0 = evolve_constant_density(
        fc, final_temperature=cooling_temperature,
        safety_factor=safety_factor)

    # evolve density and temperature according to free-fall collapse
    data = evolve_freefall(fc, final_density,
                           safety_factor=safety_factor)

    # make a plot of rho/f_H2 vs. T
    plots = pyplot.loglog(data["density"], data["temperature"],
                          color="black", label="T$_{gas}$")
    plots.extend(
        pyplot.loglog(data["density"], data["dust_temperature"],
                      color="black", linestyle="--", label="T$_{dust}$"))
    pyplot.xlabel("$\\rho$ [g/cm$^{3}$]")
    pyplot.ylabel("T [K]")

    pyplot.twinx()
    plots.extend(
        pyplot.loglog(data["density"], data["H2I"] / data["density"],
                      color="red", label="f$_{H2}$"))
    pyplot.ylabel("H$_{2}$ fraction")
    pyplot.legend(plots, [plot.get_label() for plot in plots],
                  loc="lower right")

    if os.environ.get("METAL_COOLING", 0) == "1":
        output = "freefall_metal"
    else:
        output = "freefall"

    pyplot.tight_layout()
    pyplot.savefig("%s.png" % output)

    # save data arrays as a yt dataset
    yt.save_as_dataset({}, "%s.h5" % output, data)
