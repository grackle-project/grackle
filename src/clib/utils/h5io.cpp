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

namespace {  // stuff inside an anonymous namespace is local to this file

grackle::impl::h5io::ArrayShape shape_from_space(hid_t space_id) {
  grackle::impl::h5io::ArrayShape out;
  out.ndim = -2;  // set up an output value that denotes an error
  if (space_id == H5I_INVALID_HID) {
    return out;
  }

  H5S_class_t space_type = H5Sget_simple_extent_type(space_id);
  switch (space_type) {
    case H5S_NULL: {
      out.ndim = -1;
      return out;
    }
    case H5S_SCALAR: {
      out.ndim = 0;
      return out;
    }
    case H5S_SIMPLE: {
      int ndim = H5Sget_simple_extent_ndims(space_id);
      hsize_t tmp[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION];
      if (ndim <= 0) {  // unclear what ndim of 0 would mean
        std::fprintf(stderr, "error in call to H5Sget_simple_extent_ndims\n");
        return out;
      } else if (ndim > GRACKLE_CLOUDY_TABLE_MAX_DIMENSION) {
        std::fprintf(
            stderr,
            "hdf5 array has %d dimensions, exceeding the hardcoded maximum\n",
            ndim);
        return out;
      } else if (ndim != H5Sget_simple_extent_dims(space_id, tmp, nullptr)) {
        std::fprintf(stderr, "error in call to H5Sget_simple_extent_dims\n");
        return out;
      }

      for (int i = 0; i < ndim; i++) {
        if (tmp[i] > INT64_MAX) {
          std::fprintf(stderr, "shape[%d] is too big\n", i);
          return out;
        }
        out.shape[i] = static_cast<std::int64_t>(tmp[i]);
      }
      out.ndim = ndim;  // make sure to do this last!
      return out;
    }
    default: {
      std::fprintf(stderr, "unexpected dataspace issue\n");
      return out;
    }
  }
}

}  // anonymous namespace

grackle::impl::h5io::ArrayShape grackle::impl::h5io::read_dataset_shape(
    hid_t file_id, const char* name) {
  if (name == nullptr) {
    std::fprintf(stderr, "name is a nullptr");
    return shape_from_space(H5I_INVALID_HID);
  }

  hid_t dset_id = H5Dopen(file_id, name);
  if (dset_id == H5I_INVALID_HID) {
    std::fprintf(stderr, "Failed to open dataset \"%s\".\n", name);
    return shape_from_space(H5I_INVALID_HID);
  }

  hid_t space_id = H5Dget_space(dset_id);
  H5Dclose(dset_id);

  ArrayShape out = shape_from_space(space_id);
  H5Sclose(space_id);

  if (ArrayShape_is_null(out)) {
    std::fprintf(stderr, "Error reading the dataspace of \"%s\".\n", name);
  }
  return out;
}

/// read the dataset named dset_name from file_id into buffer
int grackle::impl::h5io::read_dataset(
    hid_t file_id, const char* dset_name, double* buffer,
    const grackle::impl::h5io::ArrayShape* expected_shape) {
  if (dset_name == nullptr) {
    std::fprintf(stderr, "dset_name is a nullptr");
    return GR_FAIL;
  }

  hid_t dset_id = H5Dopen(file_id, dset_name);
  if (dset_id == H5I_INVALID_HID) {
    std::fprintf(stderr, "Failed to open dataset \"%s\".\n", dset_name);
    return GR_FAIL;
  }

  // we may perform an extra quick check about the dataset's shape
  if (expected_shape != nullptr) {
    hid_t space_id = H5Dget_space(dset_id);
    ArrayShape actual_shape = shape_from_space(space_id);
    H5Sclose(space_id);

    if (!ArrayShape_is_same(*expected_shape, actual_shape)) {
      H5Dclose(dset_id);
      std::fprintf(stderr, "The \"%s\" dataset has an unexpected shape.\n",
                   dset_name);
      return GR_FAIL;
    }
  }

  herr_t status =
      H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  H5Dclose(dset_id);
  if (status < 0) {
    std::fprintf(stderr, "Failed to read dataset \"%s\".\n", dset_name);
    return GR_FAIL;
  }

  return GR_SUCCESS;
}
