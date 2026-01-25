//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Tracks machinery to query field information
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_FIELD_INFO_DETAIL_HPP
#define GRTESTUTILS_FIELD_INFO_DETAIL_HPP

#include "grackle.h"
#include "status_reporting.h"

#include "./GrackleCtxPack.hpp"

#include <optional>
#include <string_view>

/// @defgroup FieldInfo Field Info Logic
///
/// This group all logic of the testing harness pertaining to both:
/// 1. Enumerate all fields within @ref grackle_field_data and associate
///    them with string names (to facillitate easy programatic access).
/// 2. Determine the fields Grackle expects from the current configuration.
///
/// Changes to @ref grackle_field_data should *only* directly impact the
/// test-harness logic in this group.
///
/// Longer Term Goal
/// ================
/// The longer term goal is to facillitate benchmarking with the test harness
/// logic and to make it possible to compile the test-harness logic with older
/// versions of Grackle. Certain design-decisions were made in service of this
/// goal.
///
/// Primer on Pointer-to-Data-Member
/// ================================
/// Some of this logic makes use of C++'s Pointer-to-Member functionallity (to
/// be precise, we describe a pointer to a data member). This section provides
/// a brief overview on this feature.
///
/// Just as a function-pointer is a special kind of C pointer, a
/// pointer-to-data-member is a special kind of C++ pointer.
///
/// A pointer-to-data-member is used in C++ to describe the address of a
/// struct/class's data member relative to the start of a struct. It is used
/// for to write code that in pure C might be written using `offsetof`.
///
/// For concreteness, consider the following pure C example, where `offsetof`
/// is used to access a value of the struct `S`.
/// @code{C}
/// #include <stddef.h>
/// #include <stdio.h>
/// #include <string.h>
///
/// struct S {char c, double d; double e; };
///
/// void print_e(struct S s) {
///   double tmp = 0.0;
///   memcpy(&tmp, (char*)s + offsetof(struct S, e), sizeof(double));
///   printf("s.e = %f\n", tmp);
/// }
/// @endcode
///
/// In contrast, the analogous C++ program using a pointer-to-data-member might
/// look something like the following snippet:
/// @code{C++}
/// #include <cstdio>
///
/// struct S {char c, double d; double e; };
///
/// void print_e(S s) {
///   double S::* ptr_to_member = &S::e;
///   std::printf("s.e = %f\n", s.*ptr_to_member);
/// }
/// @endcode
/// In this snippet the `ptr_to_member` variable holds a pointer-to-data-member.
/// Specifically it can hold a pointer-to-data-member that refers to any data
/// member of type `double` in the struct `S`. In other words, it could hold
/// either `&S::d` OR `&S::e`
/** @{ */

namespace grtest::field_detail {

enum class CmpOp { EQ, GEQ };

/// a requirement for enabling a field (or group of fields)
///
/// @note
/// For performance reasons, @ref param is a pointer-to-data-member. Longer
/// term, we may need to convert it to a string and use the dynamic API
struct Req {
  int chemistry_data::* param;
  CmpOp op;
  int rhs;
};

/// we implicitly assume that my_chem isn't a nullptr
inline bool check_req_(const chemistry_data& my_chem, const Req& req) {
  // ugh, this const_cast isn't great but, it's ok for our purposes
  int val = my_chem.*req.param;
  switch (req.op) {
    case CmpOp::EQ:
      return val == req.rhs;
    case CmpOp::GEQ:
      return val >= req.rhs;
    default:
      GR_INTERNAL_UNREACHABLE_ERROR();
  }
}

/// This classifies fields for the purpose of initialization
enum struct Kind { UNSET, PRIMORDIAL_SPECIES, ELEMENTWISE_SOLVER_ARG };

/// Holds information about a field member of @ref grackle_field_data
///
/// @par Further Context about Pointer-to-Data-Member
/// For the uninitiated, @ref FieldInfo::relative_addr is a special
/// kind of pointer called a pointer-to-data-member. We provide a brief primer
/// on this logic up above.
struct FieldInfo {
  gr_float* grackle_field_data::* relative_addr;  ///< address of the member
  Kind kind;                                      ///< field kind
};

/// helps implement @ref for_each_named_field
template <typename CheckReqFn, typename Fn>
constexpr void for_each_named_field_(CheckReqFn check_req, Fn fn) {
  using FInfo = FieldInfo;

#define MK_REQ(param_name, op, RHS) Req{&chemistry_data::param_name, op, RHS}

  fn("density", FInfo{&grackle_field_data::density, Kind::UNSET});
  fn("internal_energy",
     FInfo{&grackle_field_data::internal_energy, Kind::UNSET});

  if (check_req(MK_REQ(metal_cooling, CmpOp::EQ, 1))) {
    fn("metal_density", FInfo{&grackle_field_data::metal_density, Kind::UNSET});
  }

  if (check_req(MK_REQ(use_dust_density_field, CmpOp::EQ, 1))) {
    fn("dust_density", FInfo{&grackle_field_data::dust_density, Kind::UNSET});
  }

  if (check_req(MK_REQ(primordial_chemistry, CmpOp::GEQ, 1))) {
    fn("e_density",
       FInfo{&grackle_field_data::e_density, Kind::PRIMORDIAL_SPECIES});
    fn("HI_density",
       FInfo{&grackle_field_data::HI_density, Kind::PRIMORDIAL_SPECIES});
    fn("HII_density",
       FInfo{&grackle_field_data::HII_density, Kind::PRIMORDIAL_SPECIES});
    fn("HeI_density",
       FInfo{&grackle_field_data::HeI_density, Kind::PRIMORDIAL_SPECIES});
    fn("HeII_density",
       FInfo{&grackle_field_data::HeII_density, Kind::PRIMORDIAL_SPECIES});
    fn("HeIII_density",
       FInfo{&grackle_field_data::HeIII_density, Kind::PRIMORDIAL_SPECIES});
  }

  if (check_req(MK_REQ(primordial_chemistry, CmpOp::GEQ, 2))) {
    fn("HM_density",
       FInfo{&grackle_field_data::HM_density, Kind::PRIMORDIAL_SPECIES});
    fn("H2I_density",
       FInfo{&grackle_field_data::H2I_density, Kind::PRIMORDIAL_SPECIES});
    fn("H2II_density",
       FInfo{&grackle_field_data::H2II_density, Kind::PRIMORDIAL_SPECIES});
  }

  if (check_req(MK_REQ(primordial_chemistry, CmpOp::GEQ, 3))) {
    fn("DI_density",
       FInfo{&grackle_field_data::DI_density, Kind::PRIMORDIAL_SPECIES});
    fn("DII_density",
       FInfo{&grackle_field_data::DII_density, Kind::PRIMORDIAL_SPECIES});
    fn("HDI_density",
       FInfo{&grackle_field_data::HDI_density, Kind::PRIMORDIAL_SPECIES});
  }

  if (check_req(MK_REQ(primordial_chemistry, CmpOp::GEQ, 4))) {
    fn("DM_density",
       FInfo{&grackle_field_data::DM_density, Kind::PRIMORDIAL_SPECIES});
    fn("HDII_density",
       FInfo{&grackle_field_data::HDII_density, Kind::PRIMORDIAL_SPECIES});
    fn("HeHII_density",
       FInfo{&grackle_field_data::HeHII_density, Kind::PRIMORDIAL_SPECIES});
  }

  if (check_req(MK_REQ(metal_chemistry, CmpOp::EQ, 1))) {
    fn("CI_density", FInfo{&grackle_field_data::CI_density, Kind::UNSET});
    fn("CII_density", FInfo{&grackle_field_data::CII_density, Kind::UNSET});
    fn("CO_density", FInfo{&grackle_field_data::CO_density, Kind::UNSET});
    fn("CO2_density", FInfo{&grackle_field_data::CO2_density, Kind::UNSET});
    fn("OI_density", FInfo{&grackle_field_data::OI_density, Kind::UNSET});
    fn("OH_density", FInfo{&grackle_field_data::OH_density, Kind::UNSET});
    fn("H2O_density", FInfo{&grackle_field_data::H2O_density, Kind::UNSET});
    fn("O2_density", FInfo{&grackle_field_data::O2_density, Kind::UNSET});
    fn("SiI_density", FInfo{&grackle_field_data::SiI_density, Kind::UNSET});
    fn("SiOI_density", FInfo{&grackle_field_data::SiOI_density, Kind::UNSET});
    fn("SiO2I_density", FInfo{&grackle_field_data::SiO2I_density, Kind::UNSET});
    fn("CH_density", FInfo{&grackle_field_data::CH_density, Kind::UNSET});
    fn("CH2_density", FInfo{&grackle_field_data::CH2_density, Kind::UNSET});
    fn("COII_density", FInfo{&grackle_field_data::COII_density, Kind::UNSET});
    fn("OII_density", FInfo{&grackle_field_data::OII_density, Kind::UNSET});
    fn("OHII_density", FInfo{&grackle_field_data::OHII_density, Kind::UNSET});
    fn("H2OII_density", FInfo{&grackle_field_data::H2OII_density, Kind::UNSET});
    fn("H3OII_density", FInfo{&grackle_field_data::H3OII_density, Kind::UNSET});
    fn("O2II_density", FInfo{&grackle_field_data::O2II_density, Kind::UNSET});
  }

  if (check_req(MK_REQ(dust_species, CmpOp::GEQ, 1))) {
    fn("Mg_density", FInfo{&grackle_field_data::Mg_density, Kind::UNSET});
  }

  if (check_req(MK_REQ(dust_species, CmpOp::GEQ, 2))) {
    fn("Al_density", FInfo{&grackle_field_data::Al_density, Kind::UNSET});
    fn("S_density", FInfo{&grackle_field_data::S_density, Kind::UNSET});
    fn("Fe_density", FInfo{&grackle_field_data::Fe_density, Kind::UNSET});
  }

  if (check_req(MK_REQ(dust_species, CmpOp::GEQ, 1))) {
    fn("MgSiO3_dust_density",
       FInfo{&grackle_field_data::MgSiO3_dust_density, Kind::UNSET});
    fn("AC_dust_density",
       FInfo{&grackle_field_data::AC_dust_density, Kind::UNSET});
  }

  if (check_req(MK_REQ(dust_species, CmpOp::GEQ, 2))) {
    fn("SiM_dust_density",
       FInfo{&grackle_field_data::SiM_dust_density, Kind::UNSET});
    fn("FeM_dust_density",
       FInfo{&grackle_field_data::FeM_dust_density, Kind::UNSET});
    fn("Mg2SiO4_dust_density",
       FInfo{&grackle_field_data::Mg2SiO4_dust_density, Kind::UNSET});
    fn("Fe3O4_dust_density",
       FInfo{&grackle_field_data::Fe3O4_dust_density, Kind::UNSET});
    fn("SiO2_dust_density",
       FInfo{&grackle_field_data::SiO2_dust_density, Kind::UNSET});
    fn("MgO_dust_density",
       FInfo{&grackle_field_data::MgO_dust_density, Kind::UNSET});
    fn("FeS_dust_density",
       FInfo{&grackle_field_data::FeS_dust_density, Kind::UNSET});
    fn("Al2O3_dust_density",
       FInfo{&grackle_field_data::Al2O3_dust_density, Kind::UNSET});
  }

  if (check_req(MK_REQ(dust_species, CmpOp::GEQ, 3))) {
    fn("ref_org_dust_density",
       FInfo{&grackle_field_data::ref_org_dust_density, Kind::UNSET});
    fn("vol_org_dust_density",
       FInfo{&grackle_field_data::vol_org_dust_density, Kind::UNSET});
    fn("H2O_ice_dust_density",
       FInfo{&grackle_field_data::H2O_ice_dust_density, Kind::UNSET});
  }

  // we are going to ignore the case with multi_metals == 0 since I think that
  // the note in grackle_field_data is wrong in that case.
  // - That note suggests that we use one of the following fields, depending on
  //   the value of metal_abundances.
  // - In practice, I think we actually interpret the metal_density field as
  //   though it was one of the following fields.
  // - It's not worth correcting this note and handling this case since PR #480
  //   has already been proposed to eliminate the hardcoded path names
  if (check_req(MK_REQ(multi_metals, CmpOp::EQ, 1))) {
    fn("local_ISM_metal_density",
       FInfo{&grackle_field_data::local_ISM_metal_density, Kind::UNSET});
    fn("ccsn13_metal_density",
       FInfo{&grackle_field_data::ccsn13_metal_density, Kind::UNSET});
    fn("ccsn20_metal_density",
       FInfo{&grackle_field_data::ccsn20_metal_density, Kind::UNSET});
    fn("ccsn25_metal_density",
       FInfo{&grackle_field_data::ccsn25_metal_density, Kind::UNSET});
    fn("ccsn30_metal_density",
       FInfo{&grackle_field_data::ccsn30_metal_density, Kind::UNSET});
    fn("fsn13_metal_density",
       FInfo{&grackle_field_data::fsn13_metal_density, Kind::UNSET});
    fn("fsn15_metal_density",
       FInfo{&grackle_field_data::fsn15_metal_density, Kind::UNSET});
    fn("fsn50_metal_density",
       FInfo{&grackle_field_data::fsn50_metal_density, Kind::UNSET});
    fn("fsn80_metal_density",
       FInfo{&grackle_field_data::fsn80_metal_density, Kind::UNSET});
    fn("pisn170_metal_density",
       FInfo{&grackle_field_data::pisn170_metal_density, Kind::UNSET});
    fn("pisn200_metal_density",
       FInfo{&grackle_field_data::pisn200_metal_density, Kind::UNSET});
    fn("y19_metal_density",
       FInfo{&grackle_field_data::y19_metal_density, Kind::UNSET});
  }

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  if (check_req(MK_REQ(use_volumetric_heating_rate, CmpOp::EQ, 1))) {
    fn("volumetric_heating_rate",
       FInfo{&grackle_field_data::volumetric_heating_rate,
             Kind::ELEMENTWISE_SOLVER_ARG});
  }

  // specific heating rate (provide in units [erg s^-1 g^-1]
  if (check_req(MK_REQ(use_specific_heating_rate, CmpOp::EQ, 1))) {
    fn("specific_heating_rate",
       FInfo{&grackle_field_data::specific_heating_rate,
             Kind::ELEMENTWISE_SOLVER_ARG});
  }

  if (check_req(MK_REQ(use_temperature_floor, CmpOp::EQ, 1))) {
    fn("temperature_floor", FInfo{&grackle_field_data::temperature_floor,
                                  Kind::ELEMENTWISE_SOLVER_ARG});
  }

  if (check_req(MK_REQ(use_radiative_transfer, CmpOp::EQ, 1))) {
    if (check_req(MK_REQ(primordial_chemistry, CmpOp::GEQ, 1))) {
      // radiative transfer heating rate (provide in units [erg s^-1 cm^-3])
      fn("RT_heating_rate", FInfo{&grackle_field_data::RT_heating_rate,
                                  Kind::ELEMENTWISE_SOLVER_ARG});
      // radiative transfer ionization / dissociation rate fields (provided in
      // units of [1/s])
      fn("RT_HI_ionization_rate",
         FInfo{&grackle_field_data::RT_HI_ionization_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
      fn("RT_HeI_ionization_rate",
         FInfo{&grackle_field_data::RT_HeI_ionization_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
      fn("RT_HeII_ionization_rate",
         FInfo{&grackle_field_data::RT_HeII_ionization_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
    }

    if (check_req(MK_REQ(primordial_chemistry, CmpOp::GEQ, 2))) {
      fn("RT_H2_dissociation_rate",
         FInfo{&grackle_field_data::RT_H2_dissociation_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
    }

    if (check_req(MK_REQ(radiative_transfer_HDI_dissociation, CmpOp::EQ, 1))) {
      fn("RT_HDI_dissociation_rate",
         FInfo{&grackle_field_data::RT_HDI_dissociation_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
    }

    if (check_req(MK_REQ(radiative_transfer_metal_ionization, CmpOp::EQ, 1))) {
      fn("RT_CI_ionization_rate",
         FInfo{&grackle_field_data::RT_CI_ionization_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
      fn("RT_OI_ionization_rate",
         FInfo{&grackle_field_data::RT_OI_ionization_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
    }

    if (check_req(
            MK_REQ(radiative_transfer_metal_dissociation, CmpOp::EQ, 1))) {
      fn("RT_CO_dissociation_rate",
         FInfo{&grackle_field_data::RT_CO_dissociation_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
      fn("RT_OH_dissociation_rate",
         FInfo{&grackle_field_data::RT_OH_dissociation_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
      fn("RT_H2O_dissociation_rate",
         FInfo{&grackle_field_data::RT_H2O_dissociation_rate,
               Kind::ELEMENTWISE_SOLVER_ARG});
    }
  }

  // H2_self_shielding = 2
  if (check_req(MK_REQ(H2_self_shielding, CmpOp::EQ, 2))) {
    fn("H2_self_shielding_length",
       FInfo{&grackle_field_data::H2_self_shielding_length,
             Kind::ELEMENTWISE_SOLVER_ARG});
  }

  if (check_req(MK_REQ(H2_custom_shielding, CmpOp::EQ, 1))) {
    fn("H2_custom_shielding_factor",
       FInfo{&grackle_field_data::H2_custom_shielding_factor,
             Kind::ELEMENTWISE_SOLVER_ARG});
  }

  if (check_req(MK_REQ(use_isrf_field, CmpOp::EQ, 1))) {
    fn("isrf_habing",
       FInfo{&grackle_field_data::isrf_habing, Kind::ELEMENTWISE_SOLVER_ARG});
  }

#undef MK_REQ
}

/// applies the given callback function on all selected named fields
///
/// There are 2 versions of this function. These function either
/// 1. iterate over each known named field
/// 2. iterate over each named field selected by a @ref GrackleCtxPack instance
///
/// The callback function should have has a signature that resembles
/// @code{C++}
/// void my_fn(const char* field_name, const FieldInfo& f_info);
/// @endcode
/// As per usual, to maximize performance, you should pass a lambda function or
/// a callable struct rather than an ordinary function
template <typename Fn>
constexpr void for_each_named_field(Fn f) {
  return for_each_named_field_([](const Req&) -> bool { return true; }, f);
}

template <typename Fn>
void for_each_named_field(const GrackleCtxPack& ctx_pack, Fn f) {
  GR_INTERNAL_REQUIRE(ctx_pack.is_initialized(),
                      "received an uninitialized GrackleCtxPack");
  const chemistry_data* my_chemistry = ctx_pack.my_chemistry();

  auto check_req = [my_chemistry](const Req& req) -> bool {
    return check_req_(*my_chemistry, req);
  };
  return for_each_named_field_(check_req, f);
}

/// query the FieldInfo from a field name
std::optional<FieldInfo> query_field_info(std::string_view name) noexcept;

}  // namespace grtest::field_detail

/** @}*/  // end of doxygen group

#endif  // GRTESTUTILS_FIELD_INFO_DETAIL_HPP
