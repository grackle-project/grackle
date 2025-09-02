// See LICENSE file for license and copyright information

/// @file cool1d_multi_g.hpp
/// @brief Declares signature of cool1d_multi_g

// This file was initially generated automatically during conversion of the
// cool1d_multi_g function from FORTRAN to C++

#ifndef COOL1D_MULTI_G_HPP
#define COOL1D_MULTI_G_HPP

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
#include "internal_units.h"      // InternalGrUnits
#include "internal_types.hpp"    // GrainSpeciesCollection
#include "index_helper.h"        // IndexRange
// TODO: manually removed from signature (double* comp1, double* comp2)

namespace grackle::impl {

void cool1d_multi_g(int imetal, int iter,  // double* comp1, double* comp2,
                    double* edot, double* tgas, double* mmw, double* p2d,
                    double* tdust, double* metallicity, double* dust2gas,
                    double* rhoH, gr_mask_type* itmask,
                    gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
                    chemistry_data_storage* my_rates,
                    grackle_field_data* my_fields,
                    photo_rate_storage my_uvb_rates, InternalGrUnits internalu,
                    IndexRange idx_range,
                    grackle::impl::GrainSpeciesCollection grain_temperatures,
                    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
                    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
                    grackle::impl::CoolHeatScratchBuf coolingheating_buf);

};

#endif /* COOL1D_MULTI_G_HPP */
