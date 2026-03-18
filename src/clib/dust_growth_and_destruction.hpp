#ifndef DUST_GROWTH_AND_DESTRUCTION_HPP
#define DUST_GROWTH_AND_DESTRUCTION_HPP

#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "index_helper.h"
#include "internal_units.h"
#include "phys_constants.h"
#include "fortran_func_decls.h"

namespace grackle::impl {

// Calculates dust growth rates (accretion) onto grain surfaces.
// Stores the mass change dM for each cell in growth_dM array.
void dust_growth(chemistry_data* my_chemistry, grackle_field_data* my_fields,
                 InternalGrUnits internalu, IndexRange idx_range,
                 const gr_mask_type* itmask, const double* dt_value,
                 const double* t_gas,
                 double* growth_dM  // output: mass change rate for each cell
);

// Calculates dust creation rates from stellar feedback (SNe condensation).
// Stores the mass creation rate for each cell in creation_dM array.
void dust_creation(chemistry_data* my_chemistry, grackle_field_data* my_fields,
                   InternalGrUnits internalu, IndexRange idx_range,
                   const gr_mask_type* itmask, const double* dt_value,
                   double* creation_dM  // output: mass creation rate for each cell
);

// Calculates dust destruction rates from SNe shocks and thermal sputtering.
// Stores the mass change dM for each cell in destruction_dM array.
void dust_destruction(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, const double* t_gas,
    double* destruction_dM  // output: mass change rate for each cell
);

// Update the density fields using calculated mass changes.
void dust_update(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value,
    const double* growth_dM,       // input: mass change from growth
    const double* destruction_dM,  // input: mass change from destruction
    const double* creation_dM,     // input: mass change from stellar creation
    bool dryrun);

}  // namespace grackle::impl

#endif  // DUST_GROWTH_AND_DESTRUCTION_HPP