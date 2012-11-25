from grackle_defs cimport *

cdef class chemistry_data:
    cdef c_chemistry_data data
    def __cinit__(self):
        self.data = set_default_chemistry_parameters()

    property Gamma:
        def __get__(self):
            return self.data.Gamma
        def __set__(self, val):
            self.data.Gamma = val

    property use_chemistry:
        def __get__(self):
            return self.data.use_chemistry
        def __set__(self, val):
            self.data.use_chemistry = val

    property primordial_chemistry:
        def __get__(self):
            return self.data.primordial_chemistry
        def __set__(self, val):
            self.data.primordial_chemistry = val

    property metal_cooling:
        def __get__(self):
            return self.data.metal_cooling
        def __set__(self, val):
            self.data.metal_cooling = val

    property h2_on_dust:
        def __get__(self):
            return self.data.h2_on_dust
        def __set__(self, val):
            self.data.h2_on_dust = val

    property cmb_temperature_floor:
        def __get__(self):
            return self.data.cmb_temperature_floor
        def __set__(self, val):
            self.data.cmb_temperature_floor = val

    property include_metal_heating:
        def __get__(self):
            return self.data.include_metal_heating
        def __set__(self, val):
            self.data.include_metal_heating = val

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

    property photoelectric_heating:
        def __get__(self):
            return self.data.photoelectric_heating
        def __set__(self, val):
            self.data.photoelectric_heating = val

    property cloudy_table_file:
        def __get__(self):
            return self.data.cloudy_table_file
        def __set__(self, val):
            raise NotImplementedError
