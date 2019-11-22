########################################################################
#
# Python Grackle wrapper
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this 
# software.
########################################################################

from pygrackle.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    mass_hydrogen_cgs

from grackle_defs cimport *
import numpy as np
cimport numpy as np

cdef class chemistry_data:
    cdef c_chemistry_data data
    cdef c_chemistry_data_storage rates
    cdef c_code_units units

    def __cinit__(self):
        self.data = _set_default_chemistry_parameters()

    def initialize(self):
        ret =  _initialize_chemistry_data(&self.data, &self.rates, &self.units)
        if ret is None:
            raise RuntimeError("Error initializing chemistry")
        return ret

    property use_grackle:
        def __get__(self):
            return self.data.use_grackle
        def __set__(self, val):
            self.data.use_grackle = val

    property with_radiative_cooling:
        def __get__(self):
            return self.data.with_radiative_cooling
        def __set__(self, val):
            self.data.with_radiative_cooling = val

    property primordial_chemistry:
        def __get__(self):
            return self.data.primordial_chemistry
        def __set__(self, val):
            self.data.primordial_chemistry = val

    property dust_chemistry:
        def __get__(self):
            return self.data.dust_chemistry
        def __set__(self, val):
            self.data.dust_chemistry = val

    property metal_cooling:
        def __get__(self):
            return self.data.metal_cooling
        def __set__(self, val):
            self.data.metal_cooling = val

    property UVbackground:
        def __get__(self):
            return self.data.UVbackground
        def __set__(self, val):
            self.data.UVbackground = val

    property grackle_data_file:
        def __get__(self):
            return self.data.grackle_data_file
        def __set__(self, val):
            if isinstance(val, str):
                val = val.encode('utf-8')
            self.data.grackle_data_file = val

    property cmb_temperature_floor:
        def __get__(self):
            return self.data.cmb_temperature_floor
        def __set__(self, val):
            self.data.cmb_temperature_floor = val

    property Gamma:
        def __get__(self):
            return self.data.Gamma
        def __set__(self, val):
            self.data.Gamma = val

    property h2_on_dust:
        def __get__(self):
            return self.data.h2_on_dust
        def __set__(self, val):
            self.data.h2_on_dust = val

    property use_dust_density_field:
        def __get__(self):
            return self.data.use_dust_density_field
        def __set__(self, val):
            self.data.use_dust_density_field = val

    property photoelectric_heating:
        def __get__(self):
            return self.data.photoelectric_heating
        def __set__(self, val):
            self.data.photoelectric_heating = val

    property photoelectric_heating_rate:
        def __get__(self):
            return self.data.photoelectric_heating_rate
        def __set__(self, val):
            self.data.photoelectric_heating_rate = val

    property use_isrf_field:
        def __get__(self):
            return self.data.use_isrf_field
        def __set__(self, val):
            self.data.use_isrf_field = val

    property interstellar_radiation_field:
        def __get__(self):
            return self.data.interstellar_radiation_field
        def __set__(self, val):
            self.data.interstellar_radiation_field = val

    property use_volumetric_heating_rate:
        def __get__(self):
            return self.data.use_volumetric_heating_rate
        def __set__(self, val):
            self.data.use_volumetric_heating_rate = val

    property use_specific_heating_rate:
        def __get__(self):
            return self.data.use_specific_heating_rate
        def __set__(self, val):
            self.data.use_specific_heating_rate = val

    property three_body_rate:
        def __get__(self):
            return self.data.three_body_rate
        def __set__(self, val):
            self.data.three_body_rate = val

    property cie_cooling:
        def __get__(self):
            return self.data.cie_cooling
        def __set__(self, val):
            self.data.cie_cooling = val

    property h2_optical_depth_approximation:
        def __get__(self):
            return self.data.h2_optical_depth_approximation
        def __set__(self, val):
            self.data.h2_optical_depth_approximation = val

    property ih2co:
        def __get__(self):
            return self.data.ih2co
        def __set__(self, val):
            self.data.ih2co = val

    property ipiht:
        def __get__(self):
            return self.data.ipiht
        def __set__(self, val):
            self.data.ipiht = val

    property HydrogenFractionByMass:
        def __get__(self):
            return self.data.HydrogenFractionByMass
        def __set__(self, val):
            self.data.HydrogenFractionByMass = val

    property DeuteriumToHydrogenRatio:
        def __get__(self):
            return self.data.DeuteriumToHydrogenRatio
        def __set__(self, val):
            self.data.DeuteriumToHydrogenRatio = val

    property SolarMetalFractionByMass:
        def __get__(self):
            return self.data.SolarMetalFractionByMass
        def __set__(self, val):
            self.data.SolarMetalFractionByMass = val

    property local_dust_to_gas_ratio:
        def __get__(self):
            return self.data.local_dust_to_gas_ratio
        def __set__(self, val):
            self.data.local_dust_to_gas_ratio = val

    property NumberOfTemperatureBins:
        def __get__(self):
            return self.data.NumberOfTemperatureBins
        def __set__(self, val):
            self.data.NumberOfTemperatureBins = val

    property CaseBRecombination:
        def __get__(self):
            return self.data.CaseBRecombination
        def __set__(self, val):
            self.data.CaseBRecombination = val

    property TemperatureStart:
        def __get__(self):
            return self.data.TemperatureStart
        def __set__(self, val):
            self.data.TemperatureStart = val

    property TemperatureEnd:
        def __get__(self):
            return self.data.TemperatureEnd
        def __set__(self, val):
            self.data.TemperatureEnd = val

    property NumberOfDustTemperatureBins:
        def __get__(self):
            return self.data.NumberOfDustTemperatureBins
        def __set__(self, val):
            self.data.NumberOfDustTemperatureBins = val

    property DustTemperatureStart:
        def __get__(self):
            return self.data.DustTemperatureStart
        def __set__(self, val):
            self.data.DustTemperatureStart = val

    property DustTemperatureEnd:
        def __get__(self):
            return self.data.DustTemperatureEnd
        def __set__(self, val):
            self.data.DustTemperatureEnd = val

    property Compton_xray_heating:
        def __get__(self):
            return self.data.Compton_xray_heating
        def __set__(self, val):
            self.data.Compton_xray_heating = val

    property LWbackground_sawtooth_suppression:
        def __get__(self):
            return self.data.LWbackground_sawtooth_suppression
        def __set__(self, val):
            self.data.LWbackground_sawtooth_suppression = val

    property LWbackground_intensity:
        def __get__(self):
            return self.data.LWbackground_intensity
        def __set__(self, val):
            self.data.LWbackground_intensity = val

    property  UVbackground_intensity:
        def __get__(self):
            return self.data. UVbackground_intensity
        def __set__(self, val):
            self.data. UVbackground_intensity = val

    property  UVbackground_on:
        def __get__(self):
            return self.data. UVbackground_on
        def __set__(self, val):
            self.data. UVbackground_on = val

    property  UVbackground_off:
        def __get__(self):
            return self.data. UVbackground_off
        def __set__(self, val):
            self.data. UVbackground_off = val

    property  UVbackground_fullon:
        def __get__(self):
            return self.data. UVbackground_fullon
        def __set__(self, val):
            self.data. UVbackground_fullon = val

    property  UVbackground_drop:
        def __get__(self):
            return self.data. UVbackground_drop
        def __set__(self, val):
            self.data. UVbackground_drop = val

    property  cloudy_electron_fraction_factor:
        def __get__(self):
            return self.data. cloudy_electron_fraction_factor
        def __set__(self, val):
            self.data. cloudy_electron_fraction_factor = val

    property use_radiative_transfer:
        def __get__(self):
            return self.data.use_radiative_transfer
        def __set__(self, val):
            self.data.use_radiative_transfer = val

    property radiative_transfer_coupled_rate_solver:
        def __get__(self):
            return self.data.radiative_transfer_coupled_rate_solver
        def __set__(self, val):
            self.data.radiative_transfer_coupled_rate_solver = val

    property radiative_transfer_intermediate_step:
        def __get__(self):
            return self.data.radiative_transfer_intermediate_step
        def __set__(self, val):
            self.data.radiative_transfer_intermediate_step = val

    property radiative_transfer_hydrogen_only:
        def __get__(self):
            return self.data.radiative_transfer_hydrogen_only
        def __set__(self, val):
            self.data.radiative_transfer_hydrogen_only = val

    property self_shielding_method:
        def __get__(self):
            return self.data.self_shielding_method
        def __set__(self, val):
            self.data.self_shielding_method = val

    property H2_self_shielding:
        def __get__(self):
            return self.data.H2_self_shielding
        def __set__(self, val):
            self.data.H2_self_shielding = val

    property k24:
        def __get__(self):
            return self.rates.k24
        def __set__(self, val):
            self.rates.k24 = val

    property k25:
        def __get__(self):
            return self.rates.k25
        def __set__(self, val):
            self.rates.k25 = val

    property k26:
        def __get__(self):
            return self.rates.k26
        def __set__(self, val):
            self.rates.k26 = val

    property k27:
        def __get__(self):
            return self.rates.k27
        def __set__(self, val):
            self.rates.k27 = val

    property k28:
        def __get__(self):
            return self.rates.k28
        def __set__(self, val):
            self.rates.k28 = val

    property k29:
        def __get__(self):
            return self.rates.k29
        def __set__(self, val):
            self.rates.k29 = val

    property k30:
        def __get__(self):
             return self.rates.k30
        def __set__(self, val):
             self.rates.k30 = val

    property k31:
        def __get__(self):
             return self.rates.k31
        def __set__(self, val):
             self.rates.k31 = val

    property comoving_coordinates:
        def __get__(self):
            return self.units.comoving_coordinates
        def __set__(self, val):
            self.units.comoving_coordinates = val

    property density_units:
        def __get__(self):
            return self.units.density_units
        def __set__(self, val):
            self.units.density_units = val

    property length_units:
        def __get__(self):
            return self.units.length_units
        def __set__(self, val):
            self.units.length_units = val

    property time_units:
        def __get__(self):
            return self.units.time_units
        def __set__(self, val):
            self.units.time_units = val

    property velocity_units:
        def __get__(self):
            return self.units.velocity_units
        def __set__(self, val):
            self.units.velocity_units = val

    property a_units:
        def __get__(self):
            return self.units.a_units
        def __set__(self, val):
            self.units.a_units = val

    property a_value:
        def __get__(self):
            return self.units.a_value
        def __set__(self, val):
            self.units.a_value = val

    property temperature_units:
        def __get__(self):
            return mass_hydrogen_cgs * \
              self.velocity_units**2 / boltzmann_constant_cgs

    property cooling_units:
        def __get__(self):
            tbase1 = self.time_units
            if self.comoving_coordinates:
                xbase1 = self.length_units / \
                    (self.a_value * self.a_units)
                dbase1 = self.density_units * \
                    (self.a_value * self.a_units)**3
            else:
                xbase1 = self.length_units / \
                    self.a_units
                dbase1 = self.density_units * \
                    self.a_units**3

            coolunit = (self.a_units**5 * xbase1**2 *
                        mass_hydrogen_cgs**2) / (tbase1**3 * dbase1)
            return coolunit

    property energy_units:
        def __get__(self):
            return self.velocity_units**2

    property pressure_units:
        def __get__(self):
            return self.density_units * self.energy_units

cdef gr_float* get_field(fc, name):
    cdef np.ndarray rv = fc.get(name, None)
    if rv is None:
        return NULL
    else:
        return <gr_float *> rv.data

def solve_chemistry(fc, my_dt):
    cdef double dt_value = <double> my_dt

    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int grid_dimension
    grid_dimension = fc["density"].shape[0]
    cdef np.ndarray ref_gs, ref_ge
    ref_gs = np.zeros(3, dtype="int32")
    ref_ge = np.zeros(3, dtype="int32")
    ref_ge[0] = grid_dimension -1

    cdef c_field_data my_fields
    my_fields.grid_rank = 1
    my_fields.grid_dimension = &grid_dimension
    my_fields.grid_start = <int *> ref_gs.data
    my_fields.grid_end = <int *> ref_ge.data
    my_fields.density = get_field(fc, "density")
    my_fields.internal_energy = get_field(fc, "energy")
    my_fields.HI_density = get_field(fc, "HI")
    my_fields.HII_density = get_field(fc, "HII")
    my_fields.HM_density = get_field(fc, "HM")
    my_fields.HeI_density = get_field(fc, "HeI")
    my_fields.HeII_density = get_field(fc, "HeII")
    my_fields.HeIII_density = get_field(fc, "HeIII")
    my_fields.H2I_density = get_field(fc, "H2I")
    my_fields.H2II_density = get_field(fc, "H2II")
    my_fields.DI_density = get_field(fc, "DI")
    my_fields.DII_density = get_field(fc, "DII")
    my_fields.HDI_density = get_field(fc, "HDI")
    my_fields.e_density = get_field(fc, "de")
    my_fields.metal_density = get_field(fc, "metal")
    my_fields.dust_density = get_field(fc, "dust")
    my_fields.RT_heating_rate = get_field(fc, "RT_heating_rate")
    my_fields.volumetric_heating_rate = get_field(fc, "volumetric_heating_rate")
    my_fields.specific_heating_rate = get_field(fc, "specific_heating_rate")

    c_local_solve_chemistry(
        &my_chemistry,
        &my_rates,
        &my_units,
        &my_fields,
        my_dt)

def calculate_cooling_time(fc):
    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int grid_dimension
    grid_dimension = fc["density"].shape[0]
    cdef np.ndarray ref_gs, ref_ge
    ref_gs = np.zeros(3, dtype="int32")
    ref_ge = np.zeros(3, dtype="int32")
    ref_ge[0] = grid_dimension -1

    cdef c_field_data my_fields
    my_fields.grid_rank = 1
    my_fields.grid_dimension = &grid_dimension
    my_fields.grid_start = <int *> ref_gs.data
    my_fields.grid_end = <int *> ref_ge.data
    my_fields.density = get_field(fc, "density")
    my_fields.internal_energy = get_field(fc, "energy")
    my_fields.HI_density = get_field(fc, "HI")
    my_fields.HII_density = get_field(fc, "HII")
    my_fields.HM_density = get_field(fc, "HM")
    my_fields.HeI_density = get_field(fc, "HeI")
    my_fields.HeII_density = get_field(fc, "HeII")
    my_fields.HeIII_density = get_field(fc, "HeIII")
    my_fields.H2I_density = get_field(fc, "H2I")
    my_fields.H2II_density = get_field(fc, "H2II")
    my_fields.DI_density = get_field(fc, "DI")
    my_fields.DII_density = get_field(fc, "DII")
    my_fields.HDI_density = get_field(fc, "HDI")
    my_fields.e_density = get_field(fc, "de")
    my_fields.metal_density = get_field(fc, "metal")
    my_fields.dust_density = get_field(fc, "dust")
    my_fields.RT_heating_rate = get_field(fc, "RT_heating_rate")
    my_fields.volumetric_heating_rate = get_field(fc, "volumetric_heating_rate")
    my_fields.specific_heating_rate = get_field(fc, "specific_heating_rate")
    cdef gr_float *cooling_time = get_field(fc, "cooling_time")

    c_local_calculate_cooling_time(
        &my_chemistry,
        &my_rates,
        &my_units,
        &my_fields,
        cooling_time)

def calculate_gamma(fc):
    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int grid_dimension
    grid_dimension = fc["density"].shape[0]
    cdef np.ndarray ref_gs, ref_ge
    ref_gs = np.zeros(3, dtype="int32")
    ref_ge = np.zeros(3, dtype="int32")
    ref_ge[0] = grid_dimension -1

    cdef c_field_data my_fields
    my_fields.grid_rank = 1
    my_fields.grid_dimension = &grid_dimension
    my_fields.grid_start = <int *> ref_gs.data
    my_fields.grid_end = <int *> ref_ge.data
    my_fields.density = get_field(fc, "density")
    my_fields.internal_energy = get_field(fc, "energy")
    my_fields.HI_density = get_field(fc, "HI")
    my_fields.HII_density = get_field(fc, "HII")
    my_fields.HM_density = get_field(fc, "HM")
    my_fields.HeI_density = get_field(fc, "HeI")
    my_fields.HeII_density = get_field(fc, "HeII")
    my_fields.HeIII_density = get_field(fc, "HeIII")
    my_fields.H2I_density = get_field(fc, "H2I")
    my_fields.H2II_density = get_field(fc, "H2II")
    my_fields.DI_density = get_field(fc, "DI")
    my_fields.DII_density = get_field(fc, "DII")
    my_fields.HDI_density = get_field(fc, "HDI")
    my_fields.e_density = get_field(fc, "de")
    my_fields.metal_density = get_field(fc, "metal")
    my_fields.dust_density = get_field(fc, "dust")
    my_fields.RT_heating_rate = get_field(fc, "RT_heating_rate")
    my_fields.volumetric_heating_rate = get_field(fc, "volumetric_heating_rate")
    my_fields.specific_heating_rate = get_field(fc, "specific_heating_rate")
    cdef gr_float *gamma = get_field(fc, "gamma")

    c_local_calculate_gamma(
        &my_chemistry,
        &my_rates,
        &my_units,
        &my_fields,
        gamma)

def calculate_pressure(fc):
    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int grid_dimension
    grid_dimension = fc["density"].shape[0]
    cdef np.ndarray ref_gs, ref_ge
    ref_gs = np.zeros(3, dtype="int32")
    ref_ge = np.zeros(3, dtype="int32")
    ref_ge[0] = grid_dimension -1

    cdef c_field_data my_fields
    my_fields.grid_rank = 1
    my_fields.grid_dimension = &grid_dimension
    my_fields.grid_start = <int *> ref_gs.data
    my_fields.grid_end = <int *> ref_ge.data
    my_fields.density = get_field(fc, "density")
    my_fields.internal_energy = get_field(fc, "energy")
    my_fields.HI_density = get_field(fc, "HI")
    my_fields.HII_density = get_field(fc, "HII")
    my_fields.HM_density = get_field(fc, "HM")
    my_fields.HeI_density = get_field(fc, "HeI")
    my_fields.HeII_density = get_field(fc, "HeII")
    my_fields.HeIII_density = get_field(fc, "HeIII")
    my_fields.H2I_density = get_field(fc, "H2I")
    my_fields.H2II_density = get_field(fc, "H2II")
    my_fields.DI_density = get_field(fc, "DI")
    my_fields.DII_density = get_field(fc, "DII")
    my_fields.HDI_density = get_field(fc, "HDI")
    my_fields.e_density = get_field(fc, "de")
    my_fields.metal_density = get_field(fc, "metal")
    my_fields.dust_density = get_field(fc, "dust")
    my_fields.RT_heating_rate = get_field(fc, "RT_heating_rate")
    my_fields.volumetric_heating_rate = get_field(fc, "volumetric_heating_rate")
    my_fields.specific_heating_rate = get_field(fc, "specific_heating_rate")
    cdef gr_float *pressure = get_field(fc, "pressure")

    c_local_calculate_pressure(
        &my_chemistry,
        &my_rates,
        &my_units,
        &my_fields,
        pressure)

def calculate_temperature(fc):
    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int grid_dimension
    grid_dimension = fc["density"].shape[0]
    cdef np.ndarray ref_gs, ref_ge
    ref_gs = np.zeros(3, dtype="int32")
    ref_ge = np.zeros(3, dtype="int32")
    ref_ge[0] = grid_dimension -1

    cdef c_field_data my_fields
    my_fields.grid_rank = 1
    my_fields.grid_dimension = &grid_dimension
    my_fields.grid_start = <int *> ref_gs.data
    my_fields.grid_end = <int *> ref_ge.data
    my_fields.density = get_field(fc, "density")
    my_fields.internal_energy = get_field(fc, "energy")
    my_fields.HI_density = get_field(fc, "HI")
    my_fields.HII_density = get_field(fc, "HII")
    my_fields.HM_density = get_field(fc, "HM")
    my_fields.HeI_density = get_field(fc, "HeI")
    my_fields.HeII_density = get_field(fc, "HeII")
    my_fields.HeIII_density = get_field(fc, "HeIII")
    my_fields.H2I_density = get_field(fc, "H2I")
    my_fields.H2II_density = get_field(fc, "H2II")
    my_fields.DI_density = get_field(fc, "DI")
    my_fields.DII_density = get_field(fc, "DII")
    my_fields.HDI_density = get_field(fc, "HDI")
    my_fields.e_density = get_field(fc, "de")
    my_fields.metal_density = get_field(fc, "metal")
    my_fields.dust_density = get_field(fc, "dust")
    my_fields.RT_heating_rate = get_field(fc, "RT_heating_rate")
    my_fields.volumetric_heating_rate = get_field(fc, "volumetric_heating_rate")
    my_fields.specific_heating_rate = get_field(fc, "specific_heating_rate")
    cdef gr_float *temperature = get_field(fc, "temperature")

    c_local_calculate_temperature(
        &my_chemistry,
        &my_rates,
        &my_units,
        &my_fields,
        temperature)

def calculate_dust_temperature(fc):
    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int grid_dimension
    grid_dimension = fc["density"].shape[0]
    cdef np.ndarray ref_gs, ref_ge
    ref_gs = np.zeros(3, dtype="int32")
    ref_ge = np.zeros(3, dtype="int32")
    ref_ge[0] = grid_dimension -1

    cdef c_field_data my_fields
    my_fields.grid_rank = 1
    my_fields.grid_dimension = &grid_dimension
    my_fields.grid_start = <int *> ref_gs.data
    my_fields.grid_end = <int *> ref_ge.data
    my_fields.density = get_field(fc, "density")
    my_fields.internal_energy = get_field(fc, "energy")
    my_fields.HI_density = get_field(fc, "HI")
    my_fields.HII_density = get_field(fc, "HII")
    my_fields.HM_density = get_field(fc, "HM")
    my_fields.HeI_density = get_field(fc, "HeI")
    my_fields.HeII_density = get_field(fc, "HeII")
    my_fields.HeIII_density = get_field(fc, "HeIII")
    my_fields.H2I_density = get_field(fc, "H2I")
    my_fields.H2II_density = get_field(fc, "H2II")
    my_fields.DI_density = get_field(fc, "DI")
    my_fields.DII_density = get_field(fc, "DII")
    my_fields.HDI_density = get_field(fc, "HDI")
    my_fields.e_density = get_field(fc, "de")
    my_fields.metal_density = get_field(fc, "metal")
    my_fields.dust_density = get_field(fc, "dust")
    my_fields.RT_heating_rate = get_field(fc, "RT_heating_rate")
    my_fields.volumetric_heating_rate = get_field(fc, "volumetric_heating_rate")
    my_fields.specific_heating_rate = get_field(fc, "specific_heating_rate")
    cdef gr_float *dust_temperature = get_field(fc, "dust_temperature")

    c_local_calculate_dust_temperature(
        &my_chemistry,
        &my_rates,
        &my_units,
        &my_fields,
        dust_temperature)
