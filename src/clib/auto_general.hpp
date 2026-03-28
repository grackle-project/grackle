//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares the general auto-generated functions
///
//===----------------------------------------------------------------------===//

#ifndef AUTO_GENERAL_H
#define AUTO_GENERAL_H

#include <stdio.h>    // <- declares FILE
#include "grackle.h"  // <- declares get_grackle_version

/// writes compilation flags to `fp`
void auto_show_flags(FILE* fp);

/// writes configuration information to `fp`
void auto_show_config(FILE* fp);

#endif /* AUTO_GENERAL_H */
