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

from pygrackle.grackle_wrapper import \
    calculate_cooling_time, \
    calculate_gamma, \
    calculate_pressure, \
    calculate_temperature, \
    solve_chemistry

from pygrackle.utilities.misc import \
    issue_deprecation_warning

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs

_base_fluids = ["density", "metal"]
_nd_fields   = ["energy",
                "x-velocity", "y-velocity", "z-velocity",
                "temperature", "pressure",
                "gamma", "cooling_time", "mu", "nH"]

_fluid_names = {}
_fluid_names[0] = _base_fluids
_fluid_names[1] = _fluid_names[0] + \
  ["HI", "HII", "HeI", "HeII", "HeIII", "de"]
_fluid_names[2] = _fluid_names[1] + \
  ["H2I", "H2II", "HM"]
_fluid_names[3] = _fluid_names[2] + \
  ["DI", "DII", "HDI"]

_rad_trans_names = ['RT_heating_rate', 'RT_HI_ionization_rate',
                    'RT_HeI_ionization_rate', 'RT_HeII_ionization_rate',
                    'RT_H2_dissociation_rate']

_extra_fields = {}
_extra_fields[2] = ["H2_self_shielding_length"]
_extra_fields[3] = _extra_fields[2] + []

try:
    xrange
except NameError:
    xrange = range

class FluidContainer(dict):
    def __init__(self, chemistry_data, n_vals, dtype="float64",
                 itype="int64"):
        super(FluidContainer, self).__init__()
        self.dtype = dtype
        self.chemistry_data = chemistry_data
        self.n_vals = n_vals
        for fluid in _fluid_names[self.chemistry_data.primordial_chemistry] + \
        _extra_fields.get(self.chemistry_data.primordial_chemistry, []) + \
        _nd_fields:
            self._setup_fluid(fluid)
        if self.chemistry_data.use_radiative_transfer:
            for fluid in _rad_trans_names:
                self._setup_fluid(fluid)

        for htype in ["specific", "volumetric"]:
            if getattr(self.chemistry_data, "use_%s_heating_rate" % htype, 0):
                self._setup_fluid("%s_heating_rate" % htype)

    def _setup_fluid(self, fluid_name):
        self[fluid_name] = np.zeros(self.n_vals, self.dtype)

    @property
    def cooling_units(self):
        warn = 'The cooling_units attribute is deprecated.\n' + \
          'For example, instead of fc.cooling_units, ' + \
          'use fc.chemistry_data.cooling_units.'
        issue_deprecation_warning(warn)
        return self.chemistry_data.cooling_units

    @property
    def density_fields(self):
        return _fluid_names[self.chemistry_data.primordial_chemistry]

    def calculate_hydrogen_number_density(self):
        my_chemistry = self.chemistry_data
        if my_chemistry.primordial_chemistry == 0:
            self["nH"] = my_chemistry.HydrogenFractionByMass * \
              self["density"] * my_chemistry.density_units / mass_hydrogen_cgs
            return
        nH = self["HI"] + self["HII"]
        if my_chemistry.primordial_chemistry > 1:
            nH += self["HM"] + self["H2I"] + self["H2II"]
        if my_chemistry.primordial_chemistry > 2:
            nH += self["HDI"] / 2.
        self["nH"] = nH * my_chemistry.density_units / mass_hydrogen_cgs

    def calculate_mean_molecular_weight(self):
        my_chemistry = self.chemistry_data

        # If energy has been set, calculate mu from the energy
        if not (self["energy"] == 0).all():
            self.calculate_temperature()
            self["mu"] = self["temperature"] / \
                (self["energy"] * (my_chemistry.Gamma - 1.) *
                self.chemistry_data.temperature_units)
            return
    
            
        # Default to mu=1
        self["mu"] = np.ones(self["energy"].size)

        if self.chemistry_data.primordial_chemistry == 0:
            return # mu=1
        
        # Check that (chemistry) density fields have been set.
        # Allow metals to be 0
        for field in self.density_fields:
            if field == 'metal':
                continue
            if (self[field] == 0).all():
                return

        # Calculate mu from the species densities; ignore deuterium
        nden = self["metal"]/16.
        nden += self["HI"]+self["HII"]+self["de"] + \
            (self["HeI"]+self["HeII"]+self["HeIII"])/4.
            
        if self.chemistry_data.primordial_chemistry > 1:
            nden += self["HM"]+(self["H2I"]+self["H2II"])/2.
            
        self["mu"] = self["density"]/nden

        
    def calculate_cooling_time(self):
        calculate_cooling_time(self)

    def calculate_gamma(self):
        calculate_gamma(self)

    def calculate_pressure(self):
        calculate_pressure(self)

    def calculate_temperature(self):
        calculate_temperature(self)

    def solve_chemistry(self, dt):
        solve_chemistry(self, dt)

class FieldNotFound(Exception):
    def __init__(self, field):
        self.field = field

    def __str__(self):
        return "Field '%s' not found!" % ((self.field, ))

class NotAGrid(Exception):
    def __str__(self):
        return "This routine needs a yt grid."


_grackle_to_yt = {
    'density': ('gas', 'density'),
    'HI': ('gas', 'H_p0_density'),
    'HII': ('gas', 'H_p1_density'),
    'HM': ('gas', 'HM_density'),
    'HeI': ('gas', 'He_p0_density'),
    'HeII': ('gas', 'He_p1_density'),
    'HeIII': ('gas', 'He_p2_density'),
    'H2I': ('gas', 'H2_p0_density'),
    'H2II': ('gas', 'H2_p1_density'),
    'DI': ('gas', 'D_p0_density'),
    'DII': ('gas', 'D_p1_density'),
    'HDI': ('gas', 'HD_p0_density'),
    'de': ('gas', 'El_density'),
    'metal': ('gas', 'metal_density'),
    'x-velocity': ('gas', 'velocity_x'),
    'y-velocity': ('gas', 'velocity_y'),
    'z-velocity': ('gas', 'velocity_z'),
    'energy': ('gas', 'thermal_energy'),
}

_skip = ("pressure", "temperature", "cooling_time", "gamma", "mu", "nH")

_yt_to_grackle = dict((b, a) for a, b in _grackle_to_yt.items())

def _units(chemistry_data, fname):
    if fname[1].endswith("density"):
        return chemistry_data.density_units
    elif fname[1].endswith("energy"):
        energy_units = (chemistry_data.velocity_units)**2.0
        return energy_units
    elif fname[1].startswith("velocity"):
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
    ds = grid.ds
    fields = ds.derived_field_list + ds.field_list
    ni = grid.ActiveDimensions[0]
    fc = FluidContainer(chemistry_data, ni)
    for f1, f2, conv in _needed_fields(fc):
        if f2 not in fields:
            raise FieldNotFound(f2)
    for j in range(grid.ActiveDimensions[1]):
        for k in range(grid.ActiveDimensions[2]):
            for f1, f2, conv in _needed_fields(fc):
                fc[f1][:] = grid[f2][:,j,k] / conv
            yield fc
            if not update: continue
            for f1, f2, conv in _needed_fields(fc):
                grid[f2][:,j,k] = fc[f1][:] * conv
