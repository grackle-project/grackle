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

import copy
from pygrackle.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    mass_hydrogen_cgs

from libc.limits cimport INT_MAX
from .grackle_defs cimport *
import numpy as np

cdef class chemistry_data:
    cdef _wrapped_c_chemistry_data data
    cdef c_chemistry_data_storage rates
    cdef c_code_units units
    # When self.rates is uninitialized, the following is set to None. When
    # self.initialize() is called, the following is set to a deepcopy of the
    # instance of data from that moment in time. We do this because:
    # - we need the unalterred copy of self.data for deallocation of self.rates
    # - we need to let users mutate values stored in data at any time (for
    #   backwards compatibility)
    cdef object data_copy_from_init

    def __cinit__(self):
        self.data = _wrapped_c_chemistry_data()
        self.data_copy_from_init = None

    cdef void _try_uninitialize(self):
        if self.data_copy_from_init is None:
            return # nothing to uninitialize

        c_local_free_chemistry_data(
            &((<_wrapped_c_chemistry_data?>self.data_copy_from_init).data),
            &self.rates)

        self.data_copy_from_init = None

    def __del__(self):
        # only called in Python>=3.4 (in earlier versions, there will be a data
        # leak). We define this instead of __dealloc__ to ensure that
        # self.data_copy_from_init is still in a valid state
        self._try_uninitialize()

    def initialize(self):
        self._try_uninitialize() # if self.rates was already initialized,
                                 # uninitialize it to avoid a memory leak
        ret =  local_initialize_chemistry_data(&self.data.data, &self.rates,
                                               &self.units)
        if ret is None:
            raise RuntimeError("Error initializing chemistry")

        # store a deepcopy of self.data to be used to deallocate self.rates
        # later. It's important that we do this AFTER initializing self.rates
        # because self.data may be mutated within that function
        self.data_copy_from_init = copy.deepcopy(self.data)

        return ret

    def set_velocity_units(self):
        set_velocity_units(&self.units)

    def get_velocity_units(self):
        return get_velocity_units(&self.units)

    def __getattr__(self, name):
        # This has been implemented to maintain backwards compatibility.
        #
        # This method is only called whan an attribute can't be found through
        # the normal mechanism - this tries to retrieve name from self.data
        try:
            return self.data[name]
        except:
            # this method is expected to raise AttributeError when it fails
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )

    def __setattr__(self, name, value):
        # This has been implemented to maintain backwards compatibility. This
        # performs attribute assignment. It first tries to perform normal
        # attribute assignment. If that fails, try to update a value contained
        # by self.data
        #
        # Unlike with normal class objects, we can't just call
        # super(chemistry_dat, self).__setattr__(name, value) to get default
        # behavior. For extension classes, we need to manually implement the
        # default behavior of __setattr__
        #
        # in CPython, the default behavior of __setattr__ is implemented by
        # PyObject_GenericSetAttr (https://docs.python.org/3/c-api/object.html)
        #
        # For insight into how we implement this behavior, we can look to the
        # python equivalent for PyObject_GenericGetAttr that is described at
        # https://docs.python.org/3/howto/descriptor.html?highlight=descriptor%20howto#invocation-from-an-instance

        # first, check for for a data descriptor (like a property) matching
        # name in the dictionary of classes in the object’s MRO. If present, we
        # modify it instead of other attributes
        for base in type(self).__mro__:
            cls_var = vars(base).get(name, None)
            if cls_var is None:
                continue
            descr_setter = getattr(type(cls_var), '__set__', None)
            if descr_setter is not None:
                descr_setter(cls_var, self, value)
                return # early exit

        # normally, we would now check self.__dict__ for an attribute
        # called name, but that's not applicable for this class

        # finally we add our custom logic:
        try:
            self.data[name] = value
            return # early exit
        except KeyError:
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )

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

    property h2dust:
        def __get__(self):
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins*self.NumberOfDustTemperatureBins]>(<double*> self.rates.h2dust)
            return np.asarray(memview)

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

    property gas_grain:
        def __get__(self):
            if not self.dust_chemistry and not self.h2_on_dust:
                return 0
            cdef double[:] memview = <double[:self.NumberOfTemperatureBins]>(<double*> self.rates.gas_grain)
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



cdef gr_float* get_field(object fc, object name) except? NULL:
    """
    Helper function to retrieve a field's pointer from a ``FluidContainer``

    Note
    ----
    This function signature informs cython that it needs to inject code every
    time this function is called to check whether a Python Exception was raised
    when function-call returns a ``NULL`` pointer. Exceptions can arise if a
    user accidently typed ``fc[k] = vals`` instead of ``fc[k][:] = vals``.
    """
    cdef object arr = fc.get(name, None)
    if arr is None:
        return NULL

    assert arr.size == fc.n_vals # sanity check
    cdef gr_float[::1] view = arr
    return <gr_float *> &view[0]

cdef c_field_data setup_field_data(object fc, int[::1] buf,
                                   bint include_velocity) except *:
    """
    Helper function to setup an instance of the field_data struct. The
    returned struct is only valid during the lifetimes of fc and buf.
       - fc:  holds primary data used for setup
       - buf: a contiguous typed memoryview holding at least 7 elements.

    Note: The way this function is declared causes cython to generate code to
          check whether an Python Exception was raised every time this function
          gets called. The runtime overhead should be of minimal concern
    """
    # load grid size and check for error
    cdef Py_ssize_t grid_size = fc["density"].shape[0]
    if <Py_ssize_t>(INT_MAX) < grid_size:
        raise ValueError(
            ("FluidContainer stores {} values. This function is only valid "
             "for {} values or less").format(grid_size, int(INT_MAX))
        )

    # buf is used to hold memory for grid_dimension, grid_start & grid_end
    if buf.shape[0] < 7:
        raise ValueError("buf must hold at least 7 elements")

    # determine grid_dimension, grid_start and grid_end
    cdef int* grid_dimension = &buf[0]
    grid_dimension[0] = <int>(grid_size)

    buf[1:7] = 0 # initialize all remaining buffer values to zero
    cdef int* grid_start = &buf[1]
    cdef int* grid_end = &buf[4]
    grid_end[0] = <int>(grid_size) - 1

    # now initialize my_fields
    cdef c_field_data my_fields
    gr_initialize_field_data(&my_fields)
    my_fields.grid_rank = 1
    my_fields.grid_dimension = grid_dimension
    my_fields.grid_start = grid_start
    my_fields.grid_end = grid_end
    my_fields.grid_dx = -1

    my_fields.density = get_field(fc, "density")
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
    my_fields.internal_energy = get_field(fc, "energy")
    if include_velocity:
        my_fields.x_velocity = get_field(fc, "x-velocity")
        my_fields.y_velocity = get_field(fc, "y-velocity")
        my_fields.z_velocity = get_field(fc, "z-velocity")
    my_fields.volumetric_heating_rate = get_field(fc, "volumetric_heating_rate")
    my_fields.specific_heating_rate = get_field(fc, "specific_heating_rate")
    my_fields.temperature_floor = get_field(fc, "temperature_floor")
    my_fields.RT_heating_rate = get_field(fc, "RT_heating_rate")
    my_fields.RT_HI_ionization_rate = get_field(fc, "RT_HI_ionization_rate")
    my_fields.RT_HeI_ionization_rate = get_field(fc, "RT_HeI_ionization_rate")
    my_fields.RT_HeII_ionization_rate = get_field(fc, "RT_HeII_ionization_rate")
    my_fields.RT_H2_dissociation_rate = get_field(fc, "RT_H2_dissociation_rate")
    my_fields.H2_self_shielding_length = get_field(fc, "H2_self_shielding_length")
    my_fields.H2_custom_shielding_factor = get_field(fc, "H2_custom_shielding_factor")
    my_fields.isrf_habing = get_field(fc, "isrf_habing")
    return my_fields

def solve_chemistry(fc, my_dt):
    cdef double dt_value = <double> my_dt

    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int buf[7] # used for storage by members of my_fields
    cdef c_field_data my_fields = setup_field_data(fc, buf, False)

    cdef int ret = c_local_solve_chemistry(
        &my_chemistry,
        &my_rates,
        &my_units,
        &my_fields,
        my_dt)

    if (ret == GRACKLE_FAIL_VALUE):
        raise RuntimeError(f"Error occured within local_solve_chemistry")

ctypedef int (*calc_prop_fn)(c_chemistry_data*, c_chemistry_data_storage*,
			     c_code_units*, c_field_data*, gr_float*)

cdef void _calculate_helper(object fc, object field_name, object func_name,
                            calc_prop_fn func) except *:
    # handle all of the setup before calling the function:
    cdef chemistry_data chem_data = fc.chemistry_data
    cdef c_chemistry_data my_chemistry = chem_data.data.data
    cdef c_chemistry_data_storage my_rates = chem_data.rates
    cdef c_code_units my_units = chem_data.units

    cdef int buf[7] # used for storage by members of my_fields
    cdef c_field_data my_fields = setup_field_data(fc, buf, True)
    cdef gr_float * output_field = get_field(fc, field_name)

    # now actually call the function
    cdef int ret = func(&my_chemistry, &my_rates, &my_units, &my_fields,
                        output_field)
    if (ret == GRACKLE_FAIL_VALUE):
        raise RuntimeError(f"Error occured within {func_name}")
    
def calculate_cooling_time(fc):
    _calculate_helper(fc, "cooling_time", "local_calculate_cooling_time",
                      &c_local_calculate_cooling_time)

def calculate_gamma(fc):
    _calculate_helper(fc, "gamma", "local_calculate_gamma",
                      &c_local_calculate_gamma)

def calculate_pressure(fc):
    _calculate_helper(fc, "pressure", "local_calculate_pressure",
                      &c_local_calculate_pressure)

def calculate_temperature(fc):
    _calculate_helper(fc, "temperature", "local_calculate_temperature",
                      &c_local_calculate_temperature)

def calculate_dust_temperature(fc):
    _calculate_helper(fc, "dust_temperature",
                      "local_calculate_dust_temperature",
                      &c_local_calculate_dust_temperature)

def get_grackle_version():
    cdef c_grackle_version version_struct = c_get_grackle_version()
    # all members of version_struct are string literals (i.e. don't call free)
    return {"version" : version_struct.version.decode('UTF-8'),
            "branch" : version_struct.branch.decode('UTF-8'),
            "revision" : version_struct.revision.decode('UTF-8')}

cdef list _get_parameter_name_list(const char*(*get_param)(unsigned int)):
    # the arg should be param_name_int, param_name_double, or param_name_string

    cdef list out = []
    cdef const char* tmp
    cdef unsigned int i = 0
    while True:
        tmp = get_param(i)
        if tmp is NULL:
            return out
        else:
            out.append(tmp.decode('UTF-8'))
            i+=1

cdef class _wrapped_c_chemistry_data:

    cdef c_chemistry_data data
    cdef dict _string_buffers
    # Each field in data that holds a string value, may have a corresponding
    # entry in _string_buffers (the key matches the name of the field).
    #
    # - if _string_buffers DOESN'T have an entry for <field>, then data.<field>
    #   must currently store a default value. Depending on the specifics of
    #   <field>, data.<field> holds either a NULL pointer or a pointer to a
    #   string literal (stored in the C library).
    #
    # - if _string_buffers DOES have an entry for <field>, then data.<field>
    #   points to the buffer of the Python bytes object given by
    #   _string_buffers["<field>"]. In CPython, the string is always terminated
    #   by a null character (the Python interface just hides if from users)

    def __cinit__(self):
        local_initialize_chemistry_parameters(&self.data)
        self._string_buffers = {}

    def _access_struct_field(self, key, val = None):
        """
        Return a copy of the value of the field stored in the wrapped struct.

        When val is not None, this function updates the stored value first (the
        returned value is effectively a copy of val after some type coercion)
        """
        cdef c_chemistry_data* data_ptr = &self.data

        # determine whether the field needs to be updated
        update_field = (val is not None)

        # coerce key to a bytes object and then get the analogous null
        # terminated string to pass to the C functions
        coerced_key = key.encode('ascii') if isinstance(key, str) else key
        cdef char* _key = coerced_key # _key's lifetime is tied to coerced_key

        # handle the case where key is an integer
        cdef int* int_p = local_chemistry_data_access_int(data_ptr, _key)
        if int_p is not NULL:
            if update_field:
                coerced_val = int(<int?>val)
                if coerced_val != val: # prevent invalid assignment of a
                                       # non-integer floating point
                    raise TypeError("must be convertable to an int without a "
                                    "loss of information.")
                int_p[0] = <int>coerced_val
            return int(int_p[0])

        # handle the case where key is a double
        cdef double* dbl_p = local_chemistry_data_access_double(data_ptr, _key)
        if dbl_p is not NULL:
            if update_field:
                dbl_p[0] = <double?>val # <double?> cast performs type checking
            return float(dbl_p[0])

        # handle the case where key is a string
        cdef char* tmp = NULL
        cdef char** str_p = local_chemistry_data_access_string(data_ptr, _key)
        if (str_p is not NULL): # this case requires additional care

            if update_field:
                # coerce val to a bytes object & store in self._string_buffers
                if isinstance(val, str):
                    self._string_buffers[key] = bytes(val.encode('ascii'))
                elif isinstance(val, (bytes,bytearray)):
                    self._string_buffers[key] = bytes(val)
                else:
                    raise TypeError(type(val))

                assert b'\0' not in self._string_buffers[key] # sanity check
                tmp = <bytes>self._string_buffers[key]
                str_p[0] = tmp

            elif (str_p[0] is NULL):
                # handle the case where we are accessing a default field value
                # and it was initialized to a NULL pointer
                return None

            # for simplicity, just return a new copy of the contained string.
            # (this correctly handles the case where str_p[0] holds a pointer
            # to a string literal managed by the c library)
            return bytes(str_p[0])

        raise KeyError(key) # the key is not recognized

    def __getitem__(self, key):
        return self._access_struct_field(key, val = None)

    def __setitem__(self, key, value):
        if value is None:
            raise ValueError("Members of c_chemistry_data are not allowed to "
                             "be assigned values of None")
        self._access_struct_field(key, val = value)

    @staticmethod
    def int_keys(): return _get_parameter_name_list(param_name_int)

    @staticmethod
    def double_keys(): return _get_parameter_name_list(param_name_double)

    @staticmethod
    def string_keys(): return _get_parameter_name_list(param_name_string)

    def __iter__(self):
        return iter(self.int_keys() + self.double_keys() + self.string_keys())

    def __len__(self):
        return (len(self.int_keys()) + len(self.double_keys()) +
                len(self.string_keys()))

    def __deepcopy__(self, memo):
        """
        Custom deepcopy implementation. The memo dictionary is an opaque object
        """

        cdef _wrapped_c_chemistry_data out = _wrapped_c_chemistry_data()

        # the following implementation would be faster and should work, but I
        # have concerns that it could be broken easily:
        #    out.data = self.data
        #    for k in self.string_keys():
        #        out[k] = self[k]
        #    return out
        #
        # Our actual implementation is much more robust (harder to break):
        for k in self:
            out[k] = self[k]
        return out


# The following snippet exists for testing purposes. It makes use of a macro
# defined in setup.py, to specify the type of build used for libgrackle
cdef extern from *:
    """
    int traditional_in_source_build(void) {return TRADITIONAL_IN_SOURCE_BUILD;}
    """
    int traditional_in_source_build() # returns 1 for True or 0 for False

def uses_in_source_build():
    # this is only intented for testing purposes (it's subject to change)
    return bool(traditional_in_source_build())
