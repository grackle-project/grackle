########################################################################
#
# Grackle python fluid container
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this 
# software.
########################################################################

import numpy as np

_base_fluids = ["density", "metal"]
_nd_fields   = ["energy",
                "x-velocity", "y-velocity", "z-velocity",
                "temperature", "pressure",
                "gamma", "cooling_time"]
_fluid_names = {}
_fluid_names[0] = _base_fluids
_fluid_names[1] = _fluid_names[0] + \
  ["HI", "HII", "HeI", "HeII", "HeIII", "de"]
_fluid_names[2] = _fluid_names[1] + \
  ["H2I", "H2II", "HM"]
_fluid_names[3] = _fluid_names[2] + \
  ["DI", "DII", "HDI"]

class FluidContainer(dict):
    def __init__(self, chemistry_data, n_vals, dtype="float64",
                 itype="int64"):
        super(FluidContainer, self).__init__()
        self.dtype = dtype
        self.chemistry_data = chemistry_data
        self.n_vals = n_vals        
        for fluid in _fluid_names[self.chemistry_data.primordial_chemistry] + _nd_fields:
            self._setup_fluid(fluid)

    def _setup_fluid(self, fluid_name):
        self[fluid_name] = np.zeros(self.n_vals, self.dtype)

    @property
    def density_fields(self):
        return _fluid_names[self.chemistry_data.primordial_chemistry]

class FieldNotFound(Exception):
    def __init__(self, field):
        self.field = field

    def __str__(self):
        return "Field '%s' not found!" % (self.field)

class NotAGrid(Exception):
    def __str__(self):
        return "This routine needs a yt grid."

_grackle_to_yt = {
    "density": "Density",
    "HI": "HI_Density",
    "HII": "HII_Density",
    "HM": "HM_Density",
    "HeI": "HeI_Density",
    "HeII": "HeII_Density",
    "HeIII": "HeIII_Density",
    "H2I": "H2I_Density",
    "H2II": "H2II_Density",
    "DI": "DI_Density",
    "DII": "DII_Density",
    "HDI": "HDI_Density",
    "de": "Electron_Density",
    "metal": "Metal_Density",
    "x-velocity": "x-velocity",
    "y-velocity": "y-velocity",
    "z-velocity": "z-velocity",
    "energy": "ThermalEnergy",
}

_skip = ("pressure", "temperature", "cooling_time", "gamma")

_yt_to_grackle = dict((b, a) for a, b in _grackle_to_yt.items())

def _units(chemistry_data, fname):
    if fname.endswith("Density"):
        return chemistry_data.density_units
    elif fname.endswith("Energy"):
        energy_units = (chemistry_data.velocity_units)**2.0
        return energy_units
    elif fname.endswith("velocity"):
        v_units = (chemistry_data.velocity_units)
        return v_units
    else:
        raise FieldNotFound(fname)

def _needed_fields(fc):
    cd = fc.chemistry_data
    for f in sorted(fc):
        if f in _skip: continue
        f2 = _grackle_to_yt[f]
        yield f, f2, _units(cd, f2)
        
def grid_to_grackle(chemistry_data, grid, update = True):
    if not hasattr(grid, "ActiveDimensions"):
        raise RuntimeError
    pf = grid.pf
    fields = pf.h.derived_field_list + pf.h.field_list
    ni = grid.ActiveDimensions[0]
    fc = FluidContainer(chemistry_data, ni)
    for f1, f2, conv in _needed_fields(fc):
        if f2 not in fields:
            raise FieldNotFound(f2)
    for j in xrange(grid.ActiveDimensions[1]):
        for k in xrange(grid.ActiveDimensions[2]):
            for f1, f2, conv in _needed_fields(fc):
                fc[f1][:] = grid[f2][:,j,k] / conv
            yield fc
            if not update: continue
            for f1, f2, conv in _needed_fields(fc):
                grid[f2][:,j,k] = fc[f1][:] * conv
