//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Define machinery for creating GoogleTest Fixtures to help test Grackle's
/// C API.
///
//===----------------------------------------------------------------------===//

#ifndef GRTEST_FIXTURE
#define GRTEST_FIXTURE

#include <gtest/gtest.h>
// because we include gtest.h here, we should NOT include this file in any
// grtest source files (in other words, this should be a header-only file)

#include "../GrackleCtxPack.hpp"
#include "../preset.hpp"
#include "../status.hpp"

#include <utility>  // std::move

namespace grtest {

/// test-fixture that can used to run one or more tests with a chemistry data
/// configuration initialized from a given chemistry presets.
///
/// This sets up a GrackleCtxPack (where the contents are configured with
/// appropriate presets) and deallocates that memory at the end of the test.
///
/// How To Use
/// ==========
/// To make use of this fixture in a test-suite called `MyFeatureTest`, you
/// need to either:
/// 1. make a type alias (via `using` or `typedef`), named `MyFeatureTest`, of
///    the relevant instantiation of this class template, OR
/// 2. make a subclass, named `MyFeatureTest`, of the relevant instantiation of
///    this class template
template <ChemPreset chem_preset, InitialUnitPreset unit_preset>
class ConfigPresetFixture : public testing::Test {
protected:
  void SetUp() override {
    // called immediately after the constructor (but before the test-case)
    ParamConf conf = ParamConf::SimplePreset(chem_preset, unit_preset);
    std::pair<GrackleCtxPack, Status> tmp = GrackleCtxPack::create(conf);
    if (tmp.second.is_ok()) {
      pack = std::move(tmp.first);
    } else if (tmp.second.is_missing_std_file()) {
      GTEST_SKIP() << tmp.second.to_string();
    } else {
      FAIL() << tmp.second.to_string();
    }
  }

  GrackleCtxPack pack;
};

/// defines a parameterized test-fixture that can be used to run a parametrized
/// set of tests that are initialized with different chemistry data presets.
///
/// This sets up a GrackleCtxPack (where the contents are configured with
/// appropriate presets) and deallocates that memory at the end of the test.
///
/// How To Use
/// ==========
/// To make use of this fixture in a test-suite called `MyFeatureTest`, you
/// need to either:
/// 1. make a type alias (via `using` or `typedef`), named `MyFeatureTest`, of
///    this class, OR
/// 2. make a subclass, named `MyFeatureTest`, of this class
class ParametrizedConfigPresetFixture
    : public testing::TestWithParam<ParamConf> {
protected:
  void SetUp() override {
    // called immediately after the constructor (but before the test-case)
    std::pair<GrackleCtxPack, Status> tmp = GrackleCtxPack::create(GetParam());
    if (tmp.second.is_ok()) {
      pack = std::move(tmp.first);
    } else if (tmp.second.is_missing_std_file()) {
      GTEST_SKIP() << tmp.second.to_string();
    } else {
      FAIL() << tmp.second.to_string();
    }
  }

  GrackleCtxPack pack;
};

}  // namespace grtest

#endif  // GRTEST_FIXTURE
