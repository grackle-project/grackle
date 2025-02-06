// See LICENSE file for license and copyright information

/// @file internal_types.C
/// @brief Implements some types used internally by Grackle
///
/// To start out, we respect the choice to try to make this a subset of C++
/// that closely resembles C, plus (con|de)structors & namespaces
///
/// We may ultimately choose to shift to more idiomatic C++

#include <cstdlib>

#include "internal_types.hpp"
#include "utils-cpp.hpp"

// in the future, it may make sense to relocate the visitor-machinery into a
// separate file

/// represents the context for allocation
struct MemberAllocCtx_{ int nelem; };

/// implements a visitor that allocates all visited struct-members
static void visitor_allocate_member_(
  grackle::impl::MemberInfo member_info, void* member_voidp, void* ctx
) {
  int nelem = ((MemberAllocCtx_*)ctx)->nelem;

  // dispatch logic based on the struct-member's kind
  // -> As I understand it, the casting of malloc's result to the appropriate
  //    type is important to C++'s object model (which compiler's use
  //    internally to determine valid optimizations)
  switch (member_info.kind) {
    case grackle::impl::MemberKind::i64_buffer: {
      long long** member = (long long**)member_voidp;
      *member = (long long*)malloc(sizeof(long long)*nelem);
      return;
    }
    case grackle::impl::MemberKind::f64_buffer: {
      double** member = (double**)member_voidp;
      *member = (double*)malloc(sizeof(double)*nelem);
      return;
    }
    default: GRIMPL_ERROR("Encountered unexpected MemberKind");
  }
}

/// implements a visitor that frees all visited struct-members
static void visitor_cleanup_member_(
  grackle::impl::MemberInfo member_info, void* member_voidp, void* ctx
) {
  // ctx should be a nullptr

  // dispatch logic based on the struct-member's kind
  switch (member_info.kind) {
    case grackle::impl::MemberKind::i64_buffer: {
      long long** member = (long long**)member_voidp;
      GRACKLE_FREE(*member);
      return;
    }
    case grackle::impl::MemberKind::f64_buffer: {
      double** member = (double**)member_voidp;
      GRACKLE_FREE(*member);
      return;
    }
    default: GRIMPL_ERROR("Encountered unexpected MemberKind");
  }
}

// -----------------------------------------------------------------

grackle::impl::CoolHeatScratchBuf grackle::impl::new_CoolHeatScratchBuf(
  int nelem
)
{
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  CoolHeatScratchBuf out;
  MemberAllocCtx_ ctx{nelem};
  grackle::impl::visit_member_CoolHeatScratchBuf(
    &out, &visitor_allocate_member_, (void*)(&ctx)
  );
  return out;
}

void grackle::impl::drop_CoolHeatScratchBuf(CoolHeatScratchBuf* ptr)
{
  grackle::impl::visit_member_CoolHeatScratchBuf(
    ptr, &visitor_cleanup_member_, NULL
  );
}

// -----------------------------------------------------------------

grackle::impl::Cool1DMultiScratchBuf grackle::impl::new_Cool1DMultiScratchBuf(
  int nelem
)
{
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  Cool1DMultiScratchBuf out;
  MemberAllocCtx_ ctx{nelem};
  grackle::impl::visit_member_Cool1DMultiScratchBuf(
    &out, &visitor_allocate_member_, (void*)(&ctx)
  );
  return out;
}

void grackle::impl::drop_Cool1DMultiScratchBuf(
  grackle::impl::Cool1DMultiScratchBuf* ptr
)
{
  grackle::impl::visit_member_Cool1DMultiScratchBuf(
    ptr, &visitor_cleanup_member_, NULL
  );
}

// -----------------------------------------------------------------

grackle::impl::LogTLinInterpScratchBuf
grackle::impl::new_LogTLinInterpScratchBuf(int nelem)
{
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::LogTLinInterpScratchBuf out;
  MemberAllocCtx_ ctx{nelem};
  grackle::impl::visit_member_LogTLinInterpScratchBuf(
    &out, &visitor_allocate_member_, (void*)(&ctx)
  );
  return out;
}

void grackle::impl::drop_LogTLinInterpScratchBuf(
  grackle::impl::LogTLinInterpScratchBuf* ptr
)
{
  grackle::impl::visit_member_LogTLinInterpScratchBuf(
    ptr, &visitor_cleanup_member_, NULL
  );
}

// -----------------------------------------------------------------

grackle::impl::GrainSpeciesCollection grackle::impl::new_GrainSpeciesCollection(
  int nelem
)
{
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  GrainSpeciesCollection out;
  double* ptr = (double*)malloc(
    sizeof(double) * nelem * OnlyGrainSpLUT::NUM_ENTRIES
  );
  for (int i = 0; i < OnlyGrainSpLUT::NUM_ENTRIES; i++) {
    out.data[i] = ptr + (i * nelem);
  }
  return out;
}

void grackle::impl::drop_GrainSpeciesCollection(
  grackle::impl::GrainSpeciesCollection* ptr
)
{
  // since we only allocate a single pointer, we only need to call free once
  free(ptr->data[0]);
}

// -----------------------------------------------------------------

grackle::impl::SpeciesCollection grackle::impl::new_SpeciesCollection(
  int nelem
) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::SpeciesCollection out;
  double* ptr = (double*)malloc(sizeof(double) * nelem * SpLUT::NUM_ENTRIES);
  for (int i = 0; i < SpLUT::NUM_ENTRIES; i++) {
    out.data[i] = ptr + (i * nelem);
  }
  return out;
}

void grackle::impl::drop_SpeciesCollection(
  grackle::impl::SpeciesCollection *ptr
) {
  // since we only allocate a single pointer, we only need to call free once
  free(ptr->data[0]);
}

// -----------------------------------------------------------------

grackle::impl::ColRecRxnRateCollection grackle::impl::new_ColRecRxnRateCollection(
  int nelem
) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::ColRecRxnRateCollection out;
  double* ptr = (double*)malloc(
    sizeof(double) * nelem * ColRecRxnLUT::NUM_ENTRIES
  );
  for (int i = 0; i < ColRecRxnLUT::NUM_ENTRIES; i++) {
    out.data[i] = ptr + (i * nelem);
  }
  return out;
}

void grackle::impl::drop_ColRecRxnRateCollection(
  grackle::impl::ColRecRxnRateCollection *ptr
) {
  // since we only allocate a single pointer, we only need to call free once
  free(ptr->data[0]);
}

// -----------------------------------------------------------------

grackle::impl::PhotoRxnRateCollection grackle::impl::new_PhotoRxnRateCollection(
  int nelem
) {
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::PhotoRxnRateCollection out;
  MemberAllocCtx_ ctx{nelem};
  grackle::impl::visit_member_PhotoRxnRateCollection(
    &out, &visitor_allocate_member_, (void*)(&ctx)
  );
  return out;
}

void grackle::impl::drop_PhotoRxnRateCollection(
  grackle::impl::PhotoRxnRateCollection* ptr
)
{
  grackle::impl::visit_member_PhotoRxnRateCollection(
    ptr, &visitor_cleanup_member_, NULL
  );
}


// -----------------------------------------------------------------

grackle::impl::ChemHeatingRates grackle::impl::new_ChemHeatingRates(
  int nelem
)
{
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  ChemHeatingRates out;
  MemberAllocCtx_ ctx{nelem};
  grackle::impl::visit_member_ChemHeatingRates(
    &out, &visitor_allocate_member_, (void*)(&ctx)
  );
  return out;

}

void grackle::impl::drop_ChemHeatingRates(
  grackle::impl::ChemHeatingRates* ptr
)
{
  grackle::impl::visit_member_ChemHeatingRates(
    ptr, &visitor_cleanup_member_, NULL
  );
}
