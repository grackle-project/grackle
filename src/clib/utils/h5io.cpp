//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements functions to help read from hdf5 files
///
//===----------------------------------------------------------------------===//
#include <cstdio>

#include "hdf5.h"
#include "grackle.h"
#include "grackle_macros.h"

#include "h5io.hpp"
#include "../status_reporting.h"

/// read the dataset named dset_name from file_id into buffer
int grackle::impl::h5io::read_dataset(hid_t file_id, const char* dset_name,
                                      double* buffer) {
  hid_t dset_id;
  herr_t status;
  herr_t h5_error = -1;

  dset_id = H5Dopen(file_id, dset_name);
  if (dset_id == h5_error) {
    std::fprintf(stderr, "Failed to open dataset 'z'.\n");
    return GR_FAIL;
  }

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  if (status == h5_error) {
    std::fprintf(stderr, "Failed to read dataset 'z'.\n");
    return GR_FAIL;
  }

  H5Dclose(dset_id);

  return GR_SUCCESS;
}
