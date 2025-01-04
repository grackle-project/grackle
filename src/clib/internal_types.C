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

/// common function signature
typedef void modify_member_callback(double**, void*);

/// represents the context for allocation
struct MemberAllocCtx_{ int nelem; };

/// allocates the member
static void allocate_member_(double** member, void* ctx) {
  int nelem = ((MemberAllocCtx_*)ctx)->nelem;
  *member = (double*)malloc(sizeof(double)*nelem);
}

/// frees the member
static void cleanup_member_(double** member, void* ctx) {
  GRACKLE_FREE(*member);
}

// -----------------------------------------------------------------

/// Apply a function to each data member of CoolHeatScratchBuf
///
/// @note
/// If we are willing to embrace C++, then this should accept a template
/// argument (instead of a function pointer and the callback_ctx pointer)
static void for_each_coolingheating_member(
  grackle::impl::CoolHeatScratchBuf* ptr,
  modify_member_callback* fn,
  void* callback_ctx
) {
  fn(&ptr->ceHI, callback_ctx);
  fn(&ptr->ceHeI, callback_ctx);
  fn(&ptr->ceHeII, callback_ctx);
  fn(&ptr->ciHI, callback_ctx);
  fn(&ptr->ciHeI, callback_ctx);
  fn(&ptr->ciHeIS, callback_ctx);
  fn(&ptr->ciHeII, callback_ctx);
  fn(&ptr->reHII, callback_ctx);
  fn(&ptr->reHeII1, callback_ctx);
  fn(&ptr->reHeII2, callback_ctx);
  fn(&ptr->reHeIII, callback_ctx);
  fn(&ptr->brem, callback_ctx);
  fn(&ptr->cieco, callback_ctx);
  fn(&ptr->hyd01k, callback_ctx);
  fn(&ptr->h2k01, callback_ctx);
  fn(&ptr->vibh, callback_ctx);
  fn(&ptr->roth, callback_ctx);
  fn(&ptr->rotl, callback_ctx);
  fn(&ptr->gpldl, callback_ctx);
  fn(&ptr->gphdl, callback_ctx);
  fn(&ptr->hdlte, callback_ctx);
  fn(&ptr->hdlow, callback_ctx);
}

grackle::impl::CoolHeatScratchBuf grackle::impl::new_CoolHeatScratchBuf(
  int nelem
)
{
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  CoolHeatScratchBuf out;
  MemberAllocCtx_ ctx{nelem};
  for_each_coolingheating_member(&out, &allocate_member_, (void*)(&ctx));
  return out;
}

void grackle::impl::drop_CoolHeatScratchBuf(CoolHeatScratchBuf* ptr)
{
  for_each_coolingheating_member(ptr, &cleanup_member_, NULL);
}

// -----------------------------------------------------------------

grackle::impl::LogTLinInterpScratchBuf
grackle::impl::new_LogTLinInterpScratchBuf(int nelem)
{
  GRIMPL_REQUIRE(nelem > 0, "nelem must be positive");
  grackle::impl::LogTLinInterpScratchBuf out;
  out.indixe = (long long*)malloc(sizeof(long long)*nelem);
  out.t1 = (double*)malloc(sizeof(double)*nelem);
  out.t2 = (double*)malloc(sizeof(double)*nelem);
  out.logtem = (double*)malloc(sizeof(double)*nelem);
  out.tdef = (double*)malloc(sizeof(double)*nelem);
  return out;
}

void grackle::impl::drop_LogTLinInterpScratchBuf(
  grackle::impl::LogTLinInterpScratchBuf* ptr
)
{
  GRACKLE_FREE(ptr->indixe);
  GRACKLE_FREE(ptr->t1);
  GRACKLE_FREE(ptr->t2);
  GRACKLE_FREE(ptr->logtem);
  GRACKLE_FREE(ptr->tdef);
}

