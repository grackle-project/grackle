########################################################################
#
# atomic data
#
#
# Copyright (c) Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

# Solar abundances by number (as opposed to by mass) taken from the
# Cloudy photoionization code documentation (see
# https://gitlab.nublado.org/cloudy/cloudy/-/wikis/home).
#
# Note, "H" here technically refers to protium, i.e., the
# hydrogen isotope with a mass number of 1.
# When primordial_chemistry >= 3, we also follow deuterium
# (hydrogen with 1 proton and 1 neutron), which we denote
# below as "D". Deuterium does not appear in the solar
# abundance table presented by Cloudy. This is here
# mainly to allow the FluidContainer to initialize properly.
solar_abundance = {
    "H" : 1.00e+00, "He": 1.00e-01, "Li": 2.04e-09,
    "Be": 2.63e-11, "B" : 6.17e-10, "C" : 2.45e-04,
    "N" : 8.51e-05, "O" : 4.90e-04, "F" : 3.02e-08,
    "Ne": 1.00e-04, "Na": 2.14e-06, "Mg": 3.47e-05,
    "Al": 2.95e-06, "Si": 3.47e-05, "P" : 3.20e-07,
    "S" : 1.84e-05, "Cl": 1.91e-07, "Ar": 2.51e-06,
    "K" : 1.32e-07, "Ca": 2.29e-06, "Sc": 1.48e-09,
    "Ti": 1.05e-07, "V" : 1.00e-08, "Cr": 4.68e-07,
    "Mn": 2.88e-07, "Fe": 2.82e-05, "Co": 8.32e-08,
    "Ni": 1.78e-06, "Cu": 1.62e-08, "Zn": 3.98e-08,
    "D" : 0}

primordial_elements = ("H", "He", "D")
