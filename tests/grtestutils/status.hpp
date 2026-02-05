//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declare the Status class
///
//===----------------------------------------------------------------------===//
#ifndef GRTESTUTILS_STATUS_HPP
#define GRTESTUTILS_STATUS_HPP

#include <string>

namespace grtest {

namespace status_detail {

struct TmpStatus;  // forward declaration

enum class StatusKind { OK, MISSING_STD_FILE, PARAM, ADHOC };

}  // namespace status_detail

/// An instances can be used to denote context when returning from functions
///
/// @par Aside
/// Frankly, I would prefer if we used std::expected (or a simplified backport)
/// and made this only represent errors, but we can always refactor later
class Status {
  /// tracks the kind of status
  status_detail::StatusKind kind;

  /// for certain kinds of statuses, this provides extra context
  std::string s;

public:
  Status() : kind(status_detail::StatusKind::OK), s() {}

  /// this constructor is used by factory functions
  explicit Status(status_detail::TmpStatus&& tmp);

  /// checks if the status is ok
  bool is_ok() const { return kind == status_detail::StatusKind::OK; }

  /// checks if the status is an error
  bool is_err() const { return !is_ok(); }

  /// coerce to a string representation
  std::string to_string() const;

  /// Return true when `this` indicates that a standard datafile can't be found
  bool is_missing_std_file() const {
    return kind == status_detail::StatusKind::MISSING_STD_FILE;
  }
};

/// construct an Ok status
inline Status OkStatus() { return Status(); }

/// This namespace holds factory functions for @ref Status instances that
/// denote errors
namespace error {

/// Indicates that standard datafiles are missing
Status MissingStdFile();

/// There was an error setting the specified parameter name.
///
/// This means that the parameter wasn't known or the associated value had the
/// wrong type
Status Param(std::string param_name);

/// Specifies a generic string error message
Status Adhoc(std::string message);

}  // namespace error

}  // namespace grtest

#endif  // GRTESTUTILS_STATUS_HPP
