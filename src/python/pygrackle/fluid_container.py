########################################################################
#
# Grackle python fluid container
#
#
# Copyright (c) Enzo/Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np
from unyt import unyt_array

from pygrackle.grackle_wrapper import \
    calculate_cooling_time, \
    calculate_gamma, \
    calculate_pressure, \
    calculate_temperature, \
    calculate_dust_temperature, \
    solve_chemistry

from pygrackle.utilities.misc import \
    issue_deprecation_warning

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs

# approximate masses of each "element" we follow
_element_masses = {
    "H": 1,
    "D": 2,
    "He": 4,
    "e": 1,
    "metal": 16,
}

_species_info = {
    "HI": {"H": 1},
    "HII": {"H": 1},
    "HeI": {"He": 1},
    "HeII": {"He": 1},
    "HeIII": {"He": 1},
    "H2I": {"H": 2},
    "H2II": {"H": 2},
    "HM": {"H": 1},
    "DI": {"D": 1},
    "DII": {"D": 1},
    "HDI": {"H": 1, "D": 1},
    "e": {"e": 1},
    "metal": {"metal": 1},
}

_species_masses = {
    spec: sum([_element_masses[el] * num
               for el, num in info.items()])
    for spec, info in _species_info.items()
}

_deprecations = {
    "metal": "metal_density",
    "dust": "dust_density",
    "HI": "HI_density",
    "HII": "HII_density",
    "HeI": "HeI_density",
    "HeII": "HeII_density",
    "HeIII": "HeIII_density",
    "de": "e_density",
    "H2I": "H2I_density",
    "H2II": "H2II_density",
    "HM": "HM_density",
    "DI": "DI_density",
    "DII": "DII_density",
    "HDI": "HDI_density",
    "energy": "internal_energy",
    "x-velocity": "x_velocity",
    "y-velocity": "y_velocity",
    "z-velocity": "z_velocity",
    "mu": "mean_molecular_weight",
    "nH": "H_nuclei_density",
}

_base_densities = ["density"]
_base_extra_fields = \
  ["internal_energy",
   "x_velocity",
   "y_velocity",
   "z_velocity"]

# These are fields calculated using Grackle functions
# and thus arrays must be allocated for them
# Other fields calculated by the FluidContainer itself
# (for example, mean_molecular_weight) do not require this.
_calculated_fields = \
  ["cooling_time",
   "dust_temperature",
   "gamma",
   "pressure",
   "temperature"]

_primordial_chemistry_densities = {}
_primordial_chemistry_densities[0] = _base_densities
_primordial_chemistry_densities[1] = \
  _primordial_chemistry_densities[0] + \
  ["HI_density",
   "HII_density",
   "HeI_density",
   "HeII_density",
   "HeIII_density",
   "e_density"]
_primordial_chemistry_densities[2] = \
  _primordial_chemistry_densities[1] + \
  ["H2I_density",
   "H2II_density",
   "HM_density"]
_primordial_chemistry_densities[3] = \
  _primordial_chemistry_densities[2] + \
  ["DI_density",
   "DII_density",
   "HDI_density"]

_radiation_transfer_fields = \
  ["RT_heating_rate",
   "RT_HI_ionization_rate",
   "RT_HeI_ionization_rate",
   "RT_HeII_ionization_rate",
   "RT_H2_dissociation_rate"]

def _required_density_fields(my_chemistry):
    my_fields = _primordial_chemistry_densities[
        my_chemistry.primordial_chemistry].copy()
    if my_chemistry.metal_cooling == 1:
        my_fields.append("metal_density")
    if my_chemistry.dust_chemistry == 1:
        my_fields.append("dust_density")
    return my_fields

def _required_extra_fields(my_chemistry):
    my_fields = _base_extra_fields.copy()
    if my_chemistry.use_volumetric_heating_rate == 1:
        my_fields.append("volumetric_heating_rate")
    if my_chemistry.use_specific_heating_rate == 1:
        my_fields.append("specific_heating_rate")
    if my_chemistry.use_temperature_floor == 2:
        my_fields.append("temperature_floor")
    if my_chemistry.use_radiative_transfer == 1:
        my_fields.extend(_radiation_transfer_fields)
    if my_chemistry.H2_self_shielding == 2:
        my_fields.append("H2_self_shielding_length")
    if my_chemistry.H2_custom_shielding == 1:
        my_fields.append("H2_custom_shielding_factor")
    if my_chemistry.use_isrf_field == 1:
        my_fields.append("isrf_habing")
    return my_fields

class FluidContainer(dict):
    def __init__(self, chemistry_data, n_vals, dtype="float64",
                 itype="int64"):
        super(FluidContainer, self).__init__()
        self.dtype = dtype
        self.chemistry_data = chemistry_data
        self.n_vals = n_vals

        for field in self.input_fields + _calculated_fields:
            self._setup_fluid(field)

    def __getitem__(self, key):
        if key in _deprecations:
            new_field = _deprecations[key]
            warn = f"The {key} field is deprecated and will be removed in " + \
              f"Pygrackle 1.1. Use {new_field} instead."
            issue_deprecation_warning(warn)
            return self[new_field]

        return super().__getitem__(key)

    def _setup_fluid(self, fluid_name):
        self[fluid_name] = np.zeros(self.n_vals, self.dtype)

    @property
    def cooling_units(self):
        warn = "The cooling_units attribute is deprecated and will be " +\
          "removed in Pygrackle 1.1. Use chemistry_data.cooling_units instead."
        issue_deprecation_warning(warn)
        return self.chemistry_data.cooling_units

    @property
    def density_fields(self):
        return _required_density_fields(self.chemistry_data)

    @property
    def input_fields(self):
        return _required_extra_fields(self.chemistry_data) + \
          self.density_fields

    def calculate_hydrogen_number_density(self):
        warn = "calculate_hydrogen_number_density is deprecated and will " + \
          "be removed in Pygrackle 1.1. Use calculate_nuclei_density(\"H\") " + \
          "instead."
        issue_deprecation_warning(warn)
        self.calculate_nuclei_density("H")
        self["nH"] = self["H_nuclei_density"]

    def calculate_nuclei_density(self, element):
        """
        Calculate the number density of all nuclei of a given element.

        The result will be the total mass density of that element divided
        by its atomic mass. It will be stored in an entry in the
        FluidContainer named <element>_nuclei_density.

        Parameters
        ----------
        element : string
            The element in question.

        Examples
        --------
        >>> fc.calculate_nuclei_density("H")
        >>> print (fc["H_nuclei_density"])
        """
        if element not in _element_masses:
            raise ValueError(f"{element} not supported.")

        my_chemistry = self.chemistry_data
        if my_chemistry.primordial_chemistry == 0:
            rho = self["density"].copy()
            if my_chemistry.metal_cooling:
                rho -= self["metal_density"]
            if element == "H":
                rho *= my_chemistry.HydrogenFractionByMass
            elif element == "He":
                rho *= (1 - my_chemistry.HydrogenFractionByMass)
            else:
                raise ValueError(
                    f"{element} not supported for primordial_chemistry = 0.")

            self[f"{element}_nuclei_density"] = rho * my_chemistry.density_units / \
              mass_hydrogen_cgs
            return

        rho = np.zeros(self.n_vals, self.dtype)
        for spec in _species_info:
            field = f"{spec}_density"
            if field not in self:
                continue
            for my_el, my_num in _species_info[spec].items():
                if my_el != element:
                    continue
                rho += self[field] * my_num * _element_masses[my_el] / \
                  _species_masses[spec]
        self[f"{element}_nuclei_density"] = rho / _element_masses[element] * \
          my_chemistry.density_units / mass_hydrogen_cgs

    def calculate_mean_molecular_weight(self):
        """
        Calculate mean molecular weight.

        Examples
        --------
        >>> fc.calculate_mean_molecular_weight()
        >>> print (fc["mean_molecular_weight"])
        """

        # If energy has been set, calculate mu from the energy
        if not (self["internal_energy"] == 0).all():
            self.calculate_temperature()
            self.calculate_gamma()
            self["mean_molecular_weight"] = self["temperature"] / \
                (self["internal_energy"] * (self["gamma"] - 1.) *
                 self.chemistry_data.temperature_units)
            return
            
        # Default to mu=1
        self["mean_molecular_weight"] = np.ones_like(self["internal_energy"])

        if self.chemistry_data.primordial_chemistry == 0:
            return # mu=1

        n = np.zeros(self.n_vals, self.dtype)
        for field in self.density_fields:
            if field in ["density", "dust_density"]:
                continue
            spec = field[:-8]
            n += self[field] / _species_masses[spec]
        if (n == 0).any():
            print ("Warning: FluidContainer object has zero densities. "
                   "Cannot calculate a proper mean molecular weight.")
            return

        self["mean_molecular_weight"] = self["density"] / n

    def calculate_cooling_rate(self):
        """
        Calculate the cooling rate in units of erg s^-1 cm^+3.

        Examples
        --------
        >>> fc.calculate_cooling_rate()
        >>> print (fc["cooling_rate"])
        """
        self.calculate_cooling_time()

        my_chemistry = self.chemistry_data
        density_proper = self["density"] / \
            (my_chemistry.a_units *
             my_chemistry.a_value)**(3*my_chemistry.comoving_coordinates)

        cooling_rate = my_chemistry.cooling_units * self["internal_energy"] / \
          self["cooling_time"] / density_proper
        self["cooling_rate"] = cooling_rate

    def finalize_data(self):
        """
        Return field data as unyt_arrays with appropriate units.
        """

        my_chemistry = self.chemistry_data

        my_data = {}
        for field in self.density_fields:
            my_data[field] = unyt_array(
                self[field].copy() * my_chemistry.density_units, "g/cm**3")

        field_units = {
            "cooling_rate": (None, "erg*cm**3/s"),
            "cooling_time": ("time_units", "s"),
            "dust_temperature": (None, "K"),
            "internal_energy": ("energy_units", "erg/g"),
            "mean_molecular_weight": (None, ""),
            "pressure":  ("pressure_units", "dyne/cm**2"),
            "temperature": (None, "K"),
        }

        for field, (conv, units) in field_units.items():
            func = getattr(self, f"calculate_{field}", None)
            if func is not None:
                func()

            my_datum = self[field].copy()
            if conv is not None:
                my_datum *= getattr(my_chemistry, conv)
            if units:
                my_datum = unyt_array(my_datum, units)
            my_data[field] = my_datum

        return my_data

    def calculate_cooling_time(self):
        calculate_cooling_time(self)

    def calculate_gamma(self):
        calculate_gamma(self)

    def calculate_pressure(self):
        calculate_pressure(self)

    def calculate_temperature(self):
        calculate_temperature(self)

    def calculate_dust_temperature(self):
        calculate_dust_temperature(self)

    def solve_chemistry(self, dt):
        solve_chemistry(self, dt)
