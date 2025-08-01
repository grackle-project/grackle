from .grackle_defs cimport (
    c_chemistry_data, c_chemistry_data_storage, c_code_units, c_field_data
)

# see the _get_grackle_type_pack function for details
# -> there are very special rules governing the lifetime of this type (i.e. the
#    validity of its pointers)
cdef class GrackleTypePack_ExtType:
    cdef c_chemistry_data* my_chem
    cdef c_chemistry_data_storage* my_rates
    cdef c_code_units* my_units
    cdef object buffer # a numpy array used for storage by members of my_fields
    cdef c_field_data my_fields
