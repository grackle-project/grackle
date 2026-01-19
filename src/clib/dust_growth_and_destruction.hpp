#ifndef DUST_GROWTH_AND_DESTRUCTION_HPP
#define DUST_GROWTH_AND_DESTRUCTION_HPP

#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "index_helper.h"
#include "internal_units.h"
#include "phys_constants.h"
#include "fortran_func_decls.h"


namespace grackle::impl {

// Calculates and applies dust growth (accretion) onto grain surfaces.
void dust_growth(
    chemistry_data* my_chemistry,
    grackle_field_data* my_fields,
    InternalGrUnits internalu,
    IndexRange idx_range,
    const gr_mask_type* itmask,
    double* dt_value,
    double* t_gas,
    bool dryrun
);

// Calculates and applies dust destruction from SNe shocks and thermal sputtering.
void dust_destruction(
    chemistry_data* my_chemistry,
    grackle_field_data* my_fields,
    InternalGrUnits internalu,
    IndexRange idx_range,
    double* dt_value,
    double* t_gas,
    bool dryrun
);

}

#endif // DUST_GROWTH_AND_DESTRUCTION_HPP