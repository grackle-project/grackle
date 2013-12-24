########################################################################
#
# yt wrapper for grackle example script
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this 
# software.
########################################################################

from yt.mods import *
from pygrackle.api import *
from pygrackle.fluid_container import _yt_to_grackle

pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

my_chemistry = chemistry_data()
my_chemistry.use_grackle = 1
my_chemistry.primordial_chemistry = 1
my_chemistry.metal_cooling = 1
my_chemistry.grackle_data_file = "CloudyData_noUVB.h5";

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

g = pf.h.grids[0]
old = dict((f, g[f].copy()) for f in pf.h.field_list)

dt = 1e12
for fc in grid_to_grackle(my_chemistry, g):
    solve_chemistry(fc, a_value, dt)

for f in old:
    if f not in _yt_to_grackle: continue
    delta = np.abs(g[f] - old[f])/g[f]
    print f, delta.min(), delta.max()
