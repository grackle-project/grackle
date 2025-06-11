// See LICENSE file for license and copyright information

/// @file internal_types.C
/// @brief Implements some types used internally by Grackle
///
/// To start out, we respect the choice to try to make this a subset of C++
/// that closely resembles C, plus (con|de)structors & namespaces
///
/// We may ultimately choose to shift to more idiomatic C++

#include <cstdlib>  // std::malloc

#include "internal_types.hpp"
#include "utils-cpp.hpp"
#include "visitor/memory.hpp"

// -----------------------------------------------------------------

grackle::impl::CoolHeatScratchBuf grackle::impl::new_CoolHeatScratchBuf(
    int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  CoolHeatScratchBuf out;
  grackle::impl::visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  grackle::impl::visit_member(&out,
                              grackle::impl::visitor::AllocateMembers{ctx});
  return out;
}

void grackle::impl::drop_CoolHeatScratchBuf(CoolHeatScratchBuf* ptr) {
  grackle::impl::visit_member(ptr, grackle::impl::visitor::FreeMembers{});
}

// -----------------------------------------------------------------

grackle::impl::Cool1DMultiScratchBuf grackle::impl::new_Cool1DMultiScratchBuf(
    int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  Cool1DMultiScratchBuf out;
  grackle::impl::visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  grackle::impl::visit_member(&out,
                              grackle::impl::visitor::AllocateMembers{ctx});
  return out;
}

void grackle::impl::drop_Cool1DMultiScratchBuf(
    grackle::impl::Cool1DMultiScratchBuf* ptr) {
  grackle::impl::visit_member(ptr, grackle::impl::visitor::FreeMembers{});
}

// -----------------------------------------------------------------

grackle::impl::LogTLinInterpScratchBuf
grackle::impl::new_LogTLinInterpScratchBuf(int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::LogTLinInterpScratchBuf out;
  grackle::impl::visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  grackle::impl::visit_member(&out,
                              grackle::impl::visitor::AllocateMembers{ctx});
  return out;
}

void grackle::impl::drop_LogTLinInterpScratchBuf(
    grackle::impl::LogTLinInterpScratchBuf* ptr) {
  grackle::impl::visit_member(ptr, grackle::impl::visitor::FreeMembers{});
}

// -----------------------------------------------------------------

grackle::impl::GrainSpeciesCollection grackle::impl::new_GrainSpeciesCollection(
    int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  GrainSpeciesCollection out;
  double* ptr =
      (double*)malloc(sizeof(double) * nelem * OnlyGrainSpLUT::NUM_ENTRIES);
  for (int i = 0; i < OnlyGrainSpLUT::NUM_ENTRIES; i++) {
    out.data[i] = ptr + (i * nelem);
  }
  return out;
}

void grackle::impl::drop_GrainSpeciesCollection(
    grackle::impl::GrainSpeciesCollection* ptr) {
  // since we only allocate a single pointer, we only need to call free once
  free(ptr->data[0]);
}

// -----------------------------------------------------------------

grackle::impl::SpeciesCollection grackle::impl::new_SpeciesCollection(
    int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::SpeciesCollection out;
  double* ptr =
      (double*)std::malloc(sizeof(double) * nelem * SpLUT::NUM_ENTRIES);
  for (int i = 0; i < SpLUT::NUM_ENTRIES; i++) {
    out.data[i] = ptr + (i * nelem);
  }
  return out;
}

void grackle::impl::drop_SpeciesCollection(
    grackle::impl::SpeciesCollection* ptr) {
  // since we only allocate a single pointer, we only need to call free once
  free(ptr->data[0]);
}

// -----------------------------------------------------------------

grackle::impl::ColRecRxnRateCollection
grackle::impl::new_ColRecRxnRateCollection(int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::ColRecRxnRateCollection out;
  double* ptr =
      (double*)malloc(sizeof(double) * nelem * ColRecRxnLUT::NUM_ENTRIES);
  for (int i = 0; i < ColRecRxnLUT::NUM_ENTRIES; i++) {
    out.data[i] = ptr + (i * nelem);
  }
  return out;
}

void grackle::impl::drop_ColRecRxnRateCollection(
    grackle::impl::ColRecRxnRateCollection* ptr) {
  // since we only allocate a single pointer, we only need to call free once
  free(ptr->data[0]);
}

// -----------------------------------------------------------------

grackle::impl::PhotoRxnRateCollection grackle::impl::new_PhotoRxnRateCollection(
    int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::PhotoRxnRateCollection out;
  grackle::impl::visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  grackle::impl::visit_member(&out,
                              grackle::impl::visitor::AllocateMembers{ctx});
  return out;
}

void grackle::impl::drop_PhotoRxnRateCollection(
    grackle::impl::PhotoRxnRateCollection* ptr) {
  grackle::impl::visit_member(ptr, grackle::impl::visitor::FreeMembers{});
}

// -----------------------------------------------------------------

grackle::impl::ChemHeatingRates grackle::impl::new_ChemHeatingRates(int nelem) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  ChemHeatingRates out;
  grackle::impl::visitor::VisitorCtx ctx{static_cast<unsigned int>(nelem)};
  grackle::impl::visit_member(&out,
                              grackle::impl::visitor::AllocateMembers{ctx});
  return out;
}

void grackle::impl::drop_ChemHeatingRates(
    grackle::impl::ChemHeatingRates* ptr) {
  grackle::impl::visit_member(ptr, grackle::impl::visitor::FreeMembers{});
}
