########################################################################
#
# Tests rates calculated when chemistry is initialized.
#
# Written by Ewan Jones, 19/03/2021
#
########################################################################

import numpy as np 
import matplotlib.pyplot as plt 

#* Import necessary functions from grackle.
from pygrackle import chemistry_data, setup_fluid_container

#* Import necessary constants from grackle.
from pygrackle.utilities.physical_constants import mass_hydrogen_cgs

#def test_rate_initialisation():
"""
Test that the rate tables are initialized correctly.

"""
"""
#* Initialise chemistry_data instance
my_chemistry = chemistry_data()

#* Set parameters
my_chemistry.use_grackle = 1
my_chemistry.with_radiative_cooling = 0
my_chemistry.primordial_chemistry = 1
my_chemistry.metal_cooling = 0
my_chemistry.UVbackground = 0
my_chemistry.comoving_coordinates = 0
my_chemistry.a_units = 1.0
my_chemistry.a_value = 1.0
my_chemistry.density_units = mass_hydrogen_cgs
my_chemistry.length_units = 1.0
my_chemistry.time_units = 1.0
my_chemistry.velocity_units = my_chemistry.length_units / my_chemistry.time_units

#* Create fluid container from the chemistry_data instance above
my_fluidContainer = setup_fluid_container(my_chemistry, temperature=np.logspace(4.5, 9, 200),
                    converge=True, tolerance=1e-6, max_iterations=np.inf)
"""

    

            


