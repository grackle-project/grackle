//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares functions to help read from hdf5 files
///
//===----------------------------------------------------------------------===//

#ifndef UTILS_H5IO_HPP
#define UTILS_H5IO_HPP

#include "hdf5.h"

namespace grackle::impl::h5io {

/// read the dataset named dset_name from file_id into buffer
int read_dataset(hid_t file_id, const char* dset_name, double* buffer);

}  // namespace grackle::impl::h5io

#endif /* UTILS_H5IO_HPP */
