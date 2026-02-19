//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define the Status class
///
//===----------------------------------------------------------------------===//
#include "status.hpp"

#include "status_reporting.h"  // GR_INTERNAL_UNREACHABLE_ERROR
#include <utility>             // std::exchange

namespace grtest {

namespace status_detail {

/// This mirrors the internal structure of Status.
struct TmpStatus {
  StatusKind kind;
  std::string s;
};

}  // namespace status_detail

Status::Status(status_detail::TmpStatus&& tmp)
    : kind(tmp.kind), s(std::exchange(tmp.s, std::string())) {}

Status error::MissingStdFile() {
  return Status(status_detail::TmpStatus{
      status_detail::StatusKind::MISSING_STD_FILE, ""});
}

Status error::Param(std::string param_name) {
  return Status(status_detail::TmpStatus{status_detail::StatusKind::PARAM,
                                         std::move(param_name)});
}

Status error::Adhoc(std::string message) {
  return Status(status_detail::TmpStatus{status_detail::StatusKind::ADHOC,
                                         std::move(message)});
}

std::string Status::to_string() const {
  switch (kind) {
    case status_detail::StatusKind::OK:
      return "OkStatus";
    case status_detail::StatusKind::MISSING_STD_FILE: {
      return (
          "something went wrong with the test harness routine for looking up a "
          "standard datafile");
    }
    case status_detail::StatusKind::PARAM: {
      return (
          "the parameter, \"" + this->s +
          "\" isn't known to Grackle, or was associated with the wrong type");
    }
    case status_detail::StatusKind::ADHOC:
      return s;
    default:
      GR_INTERNAL_UNREACHABLE_ERROR();
  }
}

}  // namespace grtest