//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define machinery for creating GoogleTest Fixutures to help test Grackle's
/// C API.
///
//===----------------------------------------------------------------------===//

#ifndef GRTEST_FIXTURE
#define GRTEST_FIXTURE

#include <gtest/gtest.h>
// because we include gtest.h here, we should NOT include this file in any
// grtest source files (in other words, this should be a header-only file)

#include <memory>

#include "grackle.h"
#include "../preset.hpp"

namespace grtest {

/// defines a parameterized test-fixture that can used to run a parametrized
/// set of tests that are initialized with different chemistry data presets.
///
/// This sets up a GrackleCtxPack (where the contents are configured with
/// appropriate presets) and deallocates that memory at the end of the test.
///
/// How To Use
/// ==========
/// To make use of this, you might create a subclass of this type that is named
/// for the test suite. I don't love this strategy, but it seems to be the
/// standard way to do things. We can revisit this in the future.
class ParametrizedConfigPresetFixture
    : public testing::TestWithParam<FullConfPreset> {
protected:
  void SetUp() override {
    // called immediately after the constructor (but before the test-case)
    grtest::InitStatus status;
    pack_ = GrackleCtxPack::create(GetParam(), &status);
    if (!pack_.is_initialized()) {
      if (status == InitStatus::datafile_notfound) {
        GTEST_SKIP() << "something went wrong with finding the data file";
      } else {
        FAIL() << "Error in initialize_chemistry_data.";
      }
    }
  }

  GrackleCtxPack pack_;
};

}  // namespace grtest

#endif  // GRTEST_FIXTURE
