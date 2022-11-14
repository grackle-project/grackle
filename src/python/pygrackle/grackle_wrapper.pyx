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

from .grackle_defs cimport *
import numpy as np
cimport numpy as np

cdef class chemistry_data:
    cdef c_chemistry_data data
    cdef c_chemistry_data_storage rates
    cdef c_code_units units
    cdef object data_file_path

    def __cinit__(self):
        self.data = _set_default_chemistry_parameters()
        self.data_file_path = None

    def initialize(self):
        ret =  _initialize_chemistry_data(&self.data, &self.rates, &self.units)
        if ret is None:
            raise RuntimeError("Error initializing chemistry")
        return ret

    def set_velocity_units(self):
        set_velocity_units(&self.units)

    def get_velocity_units(self):
        return get_velocity_units(&self.units)

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
            # ensure that the underlying bytearray can't be modified (if it
            # grows/shrinks the `char*` allocation can be invalidated)
            return bytes(self.data.grackle_data_file)
        def __set__(self, val):
            # when Cython converts a bytearray to `char*`, the lifetime of the
            # `char*` allocation is tied to the lifetime of the original
            # bytearray object. We need to make sure that the bytearray object
            # isn't garbage collected for as long as the `char*` allocation is
            # in use. We do this by storing the bytearray as an attribute
            self.data_file_path = val
            if isinstance(self.data_file_path, str):
                self.data_file_path = self.data_file_path.encode('utf-8')
            self.data.grackle_data_file = self.data_file_path

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

    property dust_recombination_cooling:
        def __get__(self):
            return self.data.dust_recombination_cooling
        def __set__(self, val):
            self.data.dust_recombination_cooling = val

    property metal_chemistry:
        def __get__(self):
            return self.data.metal_chemistry
        def __set__(self, val):
            self.data.metal_chemistry = val

    property grain_growth:
        def __get__(self):
            return self.data.grain_growth
        def __set__(self, val):
            self.data.grain_growth = val

    property multi_metals:
        def __get__(self):
            return self.data.multi_metals
        def __set__(self, val):
            self.data.multi_metals = val

    property metal_abundances:
        def __get__(self):
            return self.data.metal_abundances
        def __set__(self, val):
            self.data.metal_abundances = val

    property dust_species:
        def __get__(self):
            return self.data.dust_species
        def __set__(self, val):
            self.data.dust_species = val

    property dust_temperature_multi:
        def __get__(self):
            return self.data.dust_temperature_multi
        def __set__(self, val):
            self.data.dust_temperature_multi = val

    property dust_sublimation:
        def __get__(self):
            return self.data.dust_sublimation
        def __set__(self, val):
            self.data.dust_sublimation = val

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

    property SN0_N:
        def __get__(self):
            return self.data.SN0_N
        def __set__(self, val):
            self.data.SN0_N = val

    property SN0_XC:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XC)

    property SN0_XO:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XO)
    
    property SN0_XMg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XMg)

    property SN0_XAl:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XAl)

    property SN0_XSi:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XSi)

    property SN0_XS:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XS)

    property SN0_XFe:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XFe)
    
    property SN0_fC:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fC)
  
    property SN0_fO:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fO)
    
    property SN0_fMg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fMg)
    
    property SN0_fAl:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fAl)
    
    property SN0_fSi:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fSi)
    
    property SN0_fS:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fS)
    
    property SN0_fFe:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fFe)
    
    property SN0_fSiM:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fSiM)
   
    property SN0_fFeM:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fFeM)
    
    property SN0_fMg2SiO4:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fMg2SiO4)
   
    property SN0_fMgSiO3:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fMgSiO3)
   
    property SN0_fFe3O4:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fFe3O4)
    
    property SN0_fAC:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fAC)
    
    property SN0_fSiO2D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fSiO2D)
    
    property SN0_fMgO:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_XC)
    
    property SN0_fFeS:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fFeS)
    
    property SN0_fAl2O3:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fAl2O3)
    
    property SN0_freforg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_freforg)
   
    property SN0_fvolorg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fvolorg)
    
    property SN0_fH2Oice:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_fH2Oice)
    
    property SN0_r0SiM:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0SiM)
    
    property SN0_r0FeM:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0FeM)
    
    property SN0_r0Mg2SiO4:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0Mg2SiO4)
    
    property SN0_r0MgSiO3:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0MgSiO3)
   
    property SN0_r0Fe3O4:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0Fe3O4)
  
    property SN0_r0AC:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0AC)
    
    property SN0_r0SiO2D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0SiO2D)
    
    property SN0_r0MgO:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0MgO)
    
    property SN0_r0FeS:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0FeS)
    
    property SN0_r0Al2O3:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0Al2O3)
    
    property SN0_r0reforg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0reforg)
    
    property SN0_r0volorg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0volorg)
    
    property SN0_r0H2Oice:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.data.SN0_r0H2Oice)    

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

    property UVbackground_redshift_on:
        def __get__(self):
            return self.data.UVbackground_redshift_on
        def __set__(self, val):
            self.data.UVbackground_redshift_on = val

    property UVbackground_redshift_off:
        def __get__(self):
            return self.data.UVbackground_redshift_off
        def __set__(self, val):
            self.data.UVbackground_redshift_off = val

    property UVbackground_redshift_fullon:
        def __get__(self):
            return self.data.UVbackground_redshift_fullon
        def __set__(self, val):
            self.data.UVbackground_redshift_fullon = val

    property UVbackground_redshift_drop:
        def __get__(self):
            return self.data.UVbackground_redshift_drop
        def __set__(self, val):
            self.data.UVbackground_redshift_drop = val

    property cloudy_electron_fraction_factor:
        def __get__(self):
            return self.data.cloudy_electron_fraction_factor
        def __set__(self, val):
            self.data.cloudy_electron_fraction_factor = val

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

    property radiative_transfer_H2II_diss:
        def __get__(self):
            return self.data.radiative_transfer_H2II_diss
        def __set__(self, val):
            self.data.radiative_transfer_H2II_diss = val

    property radiative_transfer_HDI_diss:
        def __get__(self):
            return self.data.radiative_transfer_HDI_diss
        def __set__(self, val):
            self.data.radiative_transfer_HDI_diss = val

    property radiative_transfer_metal_ion:
        def __get__(self):
            return self.data.radiative_transfer_metal_ion
        def __set__(self, val):
            self.data.radiative_transfer_metal_ion = val

    property radiative_transfer_metal_diss:
        def __get__(self):
            return self.data.radiative_transfer_metal_diss
        def __set__(self, val):
            self.data.radiative_transfer_metal_diss = val

    property radiative_transfer_use_H2_shielding:
        def __get__(self):
            return self.data.radiative_transfer_use_H2_shielding
        def __set__(self, val):
            self.data.radiative_transfer_use_H2_shielding = val

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

    property H2_custom_shielding:
        def __get__(self):
            return self.data.H2_custom_shielding
        def __set__(self, val):
            self.data.H2_custom_shielding = val

    property h2_charge_exchange_rate:
        def __get__(self):
            return self.data.h2_charge_exchange_rate
        def __set__(self, val):
            self.data.h2_charge_exchange_rate = val

    property h2_dust_rate:
        def __get__(self):
            return self.data.h2_dust_rate
        def __set__(self, val):
            self.data.h2_dust_rate = val

    property h2_h_cooling_rate:
        def __get__(self):
            return self.data.h2_h_cooling_rate
        def __set__(self, val):
            self.data.h2_h_cooling_rate = val

    property collisional_excitation_rates:
        def __get__(self):
            return self.data.collisional_excitation_rates
        def __set__(self, val):
            self.data.collisional_excitation_rates = val

    property collisional_ionisation_rates:
        def __get__(self):
            return self.data.collisional_ionisation_rates
        def __set__(self, val):
            self.data.collisional_ionisation_rates = val

    property recombination_cooling_rates:
        def __get__(self):
            return self.data.recombination_cooling_rates
        def __set__(self, val):
            self.data.recombination_cooling_rates = val

    property bremsstrahlung_cooling_rates:
        def __get__(self):
            return self.data.bremsstrahlung_cooling_rates
        def __set__(self, val):
            self.data.bremsstrahlung_cooling_rates = val

    property use_palla_salpeter_stahler_1983:
        def __get__(self):
            return self.data.use_palla_salpeter_stahler_1983
        def __set__(self, val):
            self.data.use_palla_salpeter_stahler_1983 = val

    property use_stancil_lepp_dalgarno_1998:
        def __get__(self):
            return self.data.use_stancil_lepp_dalgarno_1998
        def __set__(self, val):
            self.data.use_stancil_lepp_dalgarno_1998 = val

    property use_omukai_gas_grain:
        def __get__(self):
            return self.data.use_omukai_gas_grain
        def __set__(self, val):
            self.data.use_omukai_gas_grain = val

    property use_uniform_grain_dist_gamma_isrf:
        def __get__(self):
            return self.data.use_uniform_grain_dist_gamma_isrf
        def __set__(self, val):
            self.data.use_uniform_grain_dist_gamma_isrf = val



    property k1:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k1)
            return np.asarray(memview)

    property k2:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k2)
            return np.asarray(memview)

    property k3:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k3)
            return np.asarray(memview)

    property k4:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k4)
            return np.asarray(memview)

    property k5:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k5)
            return np.asarray(memview)

    property k6:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k6)
            return np.asarray(memview)

    property k7:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k7)
            return np.asarray(memview)
    
    property k8:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k8)
            return np.asarray(memview)

    property k9:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k9)
            return np.asarray(memview)

    property k10:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k10)
            return np.asarray(memview)

    property k11:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k11)
            return np.asarray(memview)

    property k12:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k12)
            return np.asarray(memview)

    property k13:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k13)
            return np.asarray(memview)

    property k14:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k14)
            return np.asarray(memview)

    property k15:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k15)
            return np.asarray(memview)

    property k16:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k16)
            return np.asarray(memview)

    property k17:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k17)
            return np.asarray(memview)

    property k18:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k18)
            return np.asarray(memview)

    property k19:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k19)
            return np.asarray(memview)

    property k20:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k20)
            return np.asarray(memview)

    property k21:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k21)
            return np.asarray(memview)

    property k22:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k22)
            return np.asarray(memview)

    property k23:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k23)
            return np.asarray(memview)
    
    property k13dd:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins*14]>(<double*> self.rates.k13dd)
            return np.asarray(memview)

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

    property k50:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k50)
            return np.asarray(memview)

    property k51:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k51)
            return np.asarray(memview)

    property k52:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k52)
            return np.asarray(memview)

    property k53:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k53)
            return np.asarray(memview)

    property k54:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k54)
            return np.asarray(memview)

    property k55:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k55)
            return np.asarray(memview)
    
    property k56:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k56)
            return np.asarray(memview)
        
    property k57:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k57)
            return np.asarray(memview)

    property k58:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k58)
            return np.asarray(memview)

    property k125:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k125)

    property k129:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k129)

    property k130:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k130)

    property k131:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k131)

    property k132:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k132)

    property k133:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k133)

    property k134:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k134)

    property k135:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k135)

    property k136:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k136)

    property k137:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k137)

    property k148:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k148)

    property k149:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k149)

    property k150:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k150)

    property k151:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k151)

    property k152:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k152)

    property k153:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.k153)

    property kz15:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz15)

    property kz16:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz16)

    property kz17:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz17)

    property kz18:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz18)

    property kz19:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz19)

    property kz20:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz20)

    property kz21:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz21)

    property kz22:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz22)

    property kz23:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz23)

    property kz24:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz24)

    property kz25:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz25)

    property kz26:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz26)

    property kz27:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz27)

    property kz28:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz28)

    property kz29:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz29)

    property kz30:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz30)

    property kz31:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz31)

    property kz32:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz32)

    property kz33:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz33)

    property kz34:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz34)

    property kz35:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz35)

    property kz36:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz36)

    property kz37:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz37)

    property kz38:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz38)

    property kz39:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz39)

    property kz40:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz40)

    property kz41:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz41)

    property kz42:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz42)

    property kz43:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz43)

    property kz44:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz44)

    property kz45:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz45)

    property kz46:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz46)

    property kz47:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz47)

    property kz48:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz48)

    property kz49:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz49)

    property kz50:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz50)

    property kz51:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz51)

    property kz52:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz52)

    property kz53:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz53)

    property kz54:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.kz54)

    property h2dust:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins*self.NumberOfDustTemperatureBins]>(<double*> self.rates.h2dust)
            return np.asarray(memview)

    property h2dustS:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins*self.NumberOfDustTemperatureBins]>(<double*> self.rates.h2dustS)
            return np.asarray(memview)

    property h2dustC:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins*self.NumberOfDustTemperatureBins]>(<double*> self.rates.h2dustC)
            return np.asarray(memview)

    property grain_growth_rate:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.grain_growth_rate)

    property n_cr_n:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.n_cr_n)
            return np.asarray(memview)

    property n_cr_d1:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.n_cr_d1)
            return np.asarray(memview)

    property n_cr_d2:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.n_cr_d2)
            return np.asarray(memview)

    property ceHI:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.ceHI)
            return np.asarray(memview)

    property ceHeI:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.ceHeI)
            return np.asarray(memview)

    property ceHeII:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.ceHeII)
            return np.asarray(memview)

    property ciHI:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.ciHI)
            return np.asarray(memview)

    property ciHeI:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.ciHeI)
            return np.asarray(memview)

    property ciHeIS:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.ciHeIS)
            return np.asarray(memview)
    
    property ciHeII:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.ciHeII)
            return np.asarray(memview)

    property reHII:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.reHII)
            return np.asarray(memview)

    property reHeII1:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.reHeII1)
            return np.asarray(memview)

    property reHeII2:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.reHeII2)
            return np.asarray(memview)

    property reHeIII:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.reHeIII)
            return np.asarray(memview)

    property brem:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.brem)
            return np.asarray(memview)

    property comp:
        def __get__(self):
            return self.rates.comp
        def __set__(self, val):
            self.rates.comp = val

    property comp_xray:
        def __get__(self):
            return self.rates.comp_xray
        def __set__(self, val):
            self.rates.comp_xray = val
    
    property temp_xray:
        def __get__(self):
            return self.rates.temp_xray
        def __set__(self, val):
            self.rates.temp_xray = val

    property piHI:
        def __get__(self):
            return self.rates.piHI
        def __set__(self, val):
            self.rates.piHI = val

    property piHeI:
        def __get__(self):
            return self.rates.piHeI
        def __set__(self, val):
            self.rates.piHeI = val

    property piHeII:
        def __get__(self):
            return self.rates.piHeII
        def __set__(self, val):
            self.rates.piHeII = val

    property crsHI:
        def __get__(self):
            return self.rates.crsHI
        def __set__(self, val):
            self.rates.crsHI = val

    property crsHeI:
        def __get__(self):
            return self.rates.crsHeI
        def __set__(self, val):
            self.rates.crsHeI = val

    property crsHeII:
        def __get__(self):
            return self.rates.crsHeII
        def __set__(self, val):
            self.rates.crsHeII = val

    property hyd01k:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.hyd01k)
            return np.asarray(memview)

    property h2k01:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.h2k01)
            return np.asarray(memview)
    
    property vibh:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.vibh)
            return np.asarray(memview)

    property roth:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.roth)
            return np.asarray(memview)

    property rotl:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.rotl)
            return np.asarray(memview)
    
    property GP99LowDensityLimit:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.GP99LowDensityLimit)
            return np.asarray(memview)

    property GP99HighDensityLimit:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.GP99HighDensityLimit)
            return np.asarray(memview)

    property GAHI:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.GAHI)
            return np.asarray(memview)
    
    property GAH2:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.GAH2)
            return np.asarray(memview)

    property GAHe:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.GAHe)
            return np.asarray(memview)

    property GAHp:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.GAHp)
            return np.asarray(memview)

    property GAel:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.GAel)
            return np.asarray(memview)

    property H2LTE:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.H2LTE)
            return np.asarray(memview)
    
    property HDlte:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.HDlte)
            return np.asarray(memview)

    property HDlow:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.HDlow)
            return np.asarray(memview)

    property cieco:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.cieco)
            return np.asarray(memview)

    property gammah:
        def __get__(self):
            return self.rates.gammah
        def __set__(self, val):
            self.rates.gammah = val

    property regr:
        def __get__(self):
            if not self.dust_chemistry and not self.h2_on_dust:
                return 0
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.regr)
            return np.asarray(memview)

    property gamma_isrf:
        def __get__(self):
            return self.rates.gamma_isrf
        def __set__(self, val):
            self.rates.gamma_isrf = val

    property gamma_isrf2:
        def __get__(self):
            return self.rates.gamma_isrf2
        def __set__(self, val):
            self.rates.gamma_isrf2 = val

    property gas_grain:
        def __get__(self):
            if not self.dust_chemistry and not self.h2_on_dust:
                return 0
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.gas_grain)
            return np.asarray(memview)

    property gas_grain2:
        def __get__(self):
            if not self.dust_chemistry and not self.h2_on_dust:
                return 0
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.gas_grain2)
            return np.asarray(memview)

    property cieY06:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.cieY06)
            return np.asarray(memview)

    property LH2_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LH2_N)
            return np.asarray(memview)

    property LH2_Size:
        def __get__(self):
            return self.rates.LH2_Size
        def __set__(self, val):
            self.rates.LH2_Size = val

    property LH2_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2_D)
            return np.asarray(memview)

    property LH2_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2_T)
            return np.asarray(memview)

    property LH2_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2_H)
            return np.asarray(memview)

    property LH2_dD:
        def __get__(self):
            return self.rates.LH2_dD
        def __set__(self, val):
            self.rates.LH2_dD = val

    property LH2_dT:
        def __get__(self):
            return self.rates.LH2_dT
        def __set__(self, val):
            self.rates.LH2_dT = val

    property LH2_dH:
        def __get__(self):
            return self.rates.LH2_dH
        def __set__(self, val):
            self.rates.LH2_dH = val

    property LH2_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2_L)
            return np.asarray(memview)

    property LHD_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LHD_N)
            return np.asarray(memview)

    property LHD_Size:
        def __get__(self):
            return self.rates.LHD_Size
        def __set__(self, val):
            self.rates.LHD_Size = val

    property LHD_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LHD_D)
            return np.asarray(memview)

    property LHD_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LHD_T)
            return np.asarray(memview)

    property LHD_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LHD_H)
            return np.asarray(memview)

    property LHD_dD:
        def __get__(self):
            return self.rates.LHD_dD
        def __set__(self, val):
            self.rates.LHD_dD = val

    property LHD_dT:
        def __get__(self):
            return self.rates.LHD_dT
        def __set__(self, val):
            self.rates.LHD_dT = val

    property LHD_dH:
        def __get__(self):
            return self.rates.LHD_dH
        def __set__(self, val):
            self.rates.LHD_dH = val
    
    property LHD_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LHD_L)
            return np.asarray(memview)

    property LCI_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LCI_N)
            return np.asarray(memview)

    property LCI_Size:
        def __get__(self):
            return self.rates.LCI_Size
        def __set__(self, val):
            self.rates.LCI_Size = val

    property LCI_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCI_D)
            return np.asarray(memview)

    property LCI_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCI_T)
            return np.asarray(memview)

    property LCI_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCI_H)
            return np.asarray(memview)

    property LCI_dD:
        def __get__(self):
            return self.rates.LCI_dD
        def __set__(self, val):
            self.rates.LCI_dD = val

    property LCI_dT:
        def __get__(self):
            return self.rates.LCI_dT
        def __set__(self, val):
            self.rates.LCI_dT = val

    property LCI_dH:
        def __get__(self):
            return self.rates.LCI_dH
        def __set__(self, val):
            self.rates.LCI_dH = val

    property LCI_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCI_L)
            return np.asarray(memview)

    property LCII_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LCII_N)
            return np.asarray(memview)

    property LCII_Size:
        def __get__(self):
            return self.rates.LCII_Size
        def __set__(self, val):
            self.rates.LCII_Size = val

    property LCII_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCII_D)
            return np.asarray(memview)

    property LCII_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCII_T)
            return np.asarray(memview)

    property LCII_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCII_H)
            return np.asarray(memview)

    property LCII_dD:
        def __get__(self):
            return self.rates.LCII_dD
        def __set__(self, val):
            self.rates.LCII_dD = val

    property LCII_dT:
        def __get__(self):
            return self.rates.LCII_dT
        def __set__(self, val):
            self.rates.LCII_dT = val

    property LCII_dH:
        def __get__(self):
            return self.rates.LCII_dH
        def __set__(self, val):
            self.rates.LCII_dH = val

    property LCII_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCII_L)
            return np.asarray(memview)

    property LOI_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LOI_N)
            return np.asarray(memview)

    property LOI_Size:
        def __get__(self):
            return self.rates.LOI_Size
        def __set__(self, val):
            self.rates.LOI_Size = val
  
    property LOI_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOI_D)
            return np.asarray(memview)

    property LOI_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOI_T)
            return np.asarray(memview)

    property LOI_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOI_H)
            return np.asarray(memview)

    property LOI_dD:
        def __get__(self):
            return self.rates.LOI_dD
        def __set__(self, val):
            self.rates.LOI_dD = val

    property LOI_dT:
        def __get__(self):
            return self.rates.LOI_dT
        def __set__(self, val):
            self.rates.LOI_dT = val

    property LOI_dH:
        def __get__(self):
            return self.rates.LOI_dH
        def __set__(self, val):
            self.rates.LOI_dH = val

    property LOI_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOI_L)
            return np.asarray(memview)

    property LCO_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LCO_N)
            return np.asarray(memview)

    property LCO_Size:
        def __get__(self):
            return self.rates.LCO_Size
        def __set__(self, val):
            self.rates.LCO_Size = val

    property LCO_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCO_D)
            return np.asarray(memview)

    property LCO_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCO_T)
            return np.asarray(memview)

    property LCO_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCO_H)
            return np.asarray(memview)

    property LCO_dD:
        def __get__(self):
            return self.rates.LCO_dD
        def __set__(self, val):
            self.rates.LCO_dD = val

    property LCO_dT:
        def __get__(self):
            return self.rates.LCO_dT
        def __set__(self, val):
            self.rates.LCO_dT = val

    property LCO_dH:
        def __get__(self):
            return self.rates.LCO_dH
        def __set__(self, val):
            self.rates.LCO_dH = val

    property LCO_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LCO_L)
            return np.asarray(memview)

    property LOH_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LOH_N)
            return np.asarray(memview)

    property LOH_Size:
        def __get__(self):
            return self.rates.LOH_Size
        def __set__(self, val):
            self.rates.LOH_Size = val

    property LOH_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOH_D)
            return np.asarray(memview)

    property LOH_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOH_T)
            return np.asarray(memview)

    property LOH_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOH_H)
            return np.asarray(memview)

    property LOH_dD:
        def __get__(self):
            return self.rates.LOH_dD
        def __set__(self, val):
            self.rates.LOH_dD = val

    property LOH_dT:
        def __get__(self):
            return self.rates.LOH_dT
        def __set__(self, val):
            self.rates.LOH_dT = val

    property LOH_dH:
        def __get__(self):
            return self.rates.LOH_dH
        def __set__(self, val):
            self.rates.LOH_dH = val

    property LOH_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LOH_L)
            return np.asarray(memview)

    property LH2O_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.LH2O_N)
            return np.asarray(memview)

    property LH2O_Size:
        def __get__(self):
            return self.rates.LH2O_Size
        def __set__(self, val):
            self.rates.LH2O_Size = val

    property LH2O_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2O_D)
            return np.asarray(memview)

    property LH2O_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2O_T)
            return np.asarray(memview)

    property LH2O_H:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2O_H)
            return np.asarray(memview)

    property LH2O_dD:
        def __get__(self):
            return self.rates.LH2O_dD
        def __set__(self, val):
            self.rates.LH2O_dD = val

    property LH2O_dT:
        def __get__(self):
            return self.rates.LH2O_dT
        def __set__(self, val):
            self.rates.LH2O_dT = val

    property LH2O_dH:
        def __get__(self):
            return self.rates.LH2O_dH
        def __set__(self, val):
            self.rates.LH2O_dH = val

    property LH2O_L:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.LH2O_L)
            return np.asarray(memview)

    property alphap_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.alphap_N)
            return np.asarray(memview)

    property alphap_Size:
        def __get__(self):
            return self.rates.alphap_Size
        def __set__(self, val):
            self.rates.alphap_Size = val

    property alphap_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.alphap_D)
            return np.asarray(memview)

    property alphap_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.alphap_T)
            return np.asarray(memview)

    property alphap_dD:
        def __get__(self):
            return self.rates.alphap_dD
        def __set__(self, val):
            self.rates.alphap_dD = val

    property alphap_dT:
        def __get__(self):
            return self.rates.alphap_dT
        def __set__(self, val):
            self.rates.alphap_dT = val

    property alphap_Data:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.alphap_Data)
            return np.asarray(memview)

    property grain_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.grain_N)
            return np.asarray(memview)

    property grain_Size:
        def __get__(self):
            return self.rates.grain_Size
        def __set__(self, val):
            self.rates.grain_Size = val

    property grain_D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.grain_D)
            return np.asarray(memview)

    property grain_T:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.grain_T)
            return np.asarray(memview)

    property grain_dD:
        def __get__(self):
            return self.rates.grain_dD
        def __set__(self, val):
            self.rates.grain_dD = val

    property grain_dT:
        def __get__(self):
            return self.rates.grain_dT
        def __set__(self, val):
            self.rates.grain_dT = val

    property Hgrain:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.Hgrain)
            return np.asarray(memview)

    property Tgrain:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.Tgrain)
            return np.asarray(memview)

    property Ograin:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.Ograin)
            return np.asarray(memview)

    property Lgrain:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.Lgrain)
            return np.asarray(memview)

    property gr_N:
        def __get__(self):
            cdef int[:] memview = <int[:self.NumberOfTemperatureBins]>(<int*> self.rates.gr_N)
            return np.asarray(memview)

    property gr_Size:
        def __get__(self):
            return self.rates.gr_Size
        def __set__(self, val):
            self.rates.gr_Size = val

    property gr_dT:
        def __get__(self):
            return self.rates.gr_dT
        def __set__(self, val):
            self.rates.gr_dT = val

    property gr_Td:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.gr_Td)
            return np.asarray(memview)

    property SN0_kpSiM:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpSiM)
            return np.asarray(memview)

    property SN0_kpFeM:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpFeM)
            return np.asarray(memview)

    property SN0_kpMg2SiO4:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpMg2SiO4)
            return np.asarray(memview)

    property SN0_kpMgSiO3:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpMgSiO3)
            return np.asarray(memview)

    property SN0_kpFe3O4:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpFe3O4)
            return np.asarray(memview)

    property SN0_kpAC:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpAC)
            return np.asarray(memview)

    property SN0_kpSiO2D:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpSiO2D)
            return np.asarray(memview)

    property SN0_kpMgO:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpMgO)
            return np.asarray(memview)

    property SN0_kpFeS:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpFeS)
            return np.asarray(memview)

    property SN0_kpAl2O3:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpAl2O3)
            return np.asarray(memview)

    property SN0_kpreforg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpreforg)
            return np.asarray(memview)

    property SN0_kpvolorg:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpvolorg)
            return np.asarray(memview)

    property SN0_kpH2Oice:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.SN0_kpH2Oice)
            return np.asarray(memview)

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
            return get_temperature_units(&self.units)

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
    my_fields.x_velocity = get_field(fc, "x-velocity")
    my_fields.y_velocity = get_field(fc, "y-velocity")
    my_fields.z_velocity = get_field(fc, "z-velocity")
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
    my_fields.x_velocity = get_field(fc, "x-velocity")
    my_fields.y_velocity = get_field(fc, "y-velocity")
    my_fields.z_velocity = get_field(fc, "z-velocity")
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
    my_fields.x_velocity = get_field(fc, "x-velocity")
    my_fields.y_velocity = get_field(fc, "y-velocity")
    my_fields.z_velocity = get_field(fc, "z-velocity")
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
    my_fields.x_velocity = get_field(fc, "x-velocity")
    my_fields.y_velocity = get_field(fc, "y-velocity")
    my_fields.z_velocity = get_field(fc, "z-velocity")
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
    my_fields.x_velocity = get_field(fc, "x-velocity")
    my_fields.y_velocity = get_field(fc, "y-velocity")
    my_fields.z_velocity = get_field(fc, "z-velocity")
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

def get_grackle_version():
    cdef c_grackle_version version_struct = c_get_grackle_version()
    # all members of version_struct are string literals (i.e. don't call free)
    return {"version" : version_struct.version.decode('UTF-8'),
            "branch" : version_struct.branch.decode('UTF-8'),
            "revision" : version_struct.revision.decode('UTF-8')}
