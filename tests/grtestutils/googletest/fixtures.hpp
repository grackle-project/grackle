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

// C Pre-Processor weirdness requires an additional level of indirection
// if we want B to be passed __LINE__ and actually expand to a line number
#define GRTest_CONCAT_TOKEN_HELPER_(a, b) a##b
#define GRTest_CONCAT_TOKEN_(a, b) GRTest_CONCAT_TOKEN_HELPER_(a, b)

/// Construct a GrackleCtxPack with googletest error checking
///
/// In more detail, this calls `GrackleCtxPack::create(conf)`, performs
/// associated error-checking and then stores the result in `lhs`. This is
/// roughly analogous to
/// \code{C++}
/// RsltPairT rslt_pair = ::grtest::GrackleCtxPack::create(conf);
/// if (!rslt_pair.second.is_err()) {
///   /* handle the error with Googletest machinery */
/// }
/// lhs = std::move(rslt_pair.first);
/// \endcode
/// In the above pseudocode, 2 simplifications were made:
/// - `RsltPairT` replaced std::pair<::grtest::GrackleCtxPack, ::grtest::Status>
/// - the variable name `rslt_pair` replaces a more unique variable name
///
/// @param[in, out] lhs If construction of the pack succeeds, then this will
///     be on the left hand side of the final assignment expression
/// @param[in] conf Expands to eSimpleConfPreset
///
/// @par Motivation
/// This exists so that one-off test-cases don't need to create a full-blown
/// fixture to check a given configuration (they can use this macro instead)
///
/// @warning
/// Because the error-checking logic may trigger skipping of a test or a test
/// failure, this **MUST** be called within the body of a test-case or in the
/// `SetUp()`/`TearDown()` method of a fixture. In other words, **DO NOT** call
/// this in a helper function. For more context, see:
/// https://google.github.io/googletest/advanced.html#assertion-placement
///
/// @par Broader Thoughts
/// I really don't like this. But, it currently seems like the most pragmatic
/// solution for now...
/// - Part of the issue is that I have ambitions to reuse the GrackleCtxPack
///   for benchmarking and possibly fuzz/property testing and I didn't want to
///   pigeonhole ourselves...
/// - Now that error-handling of `GrackleCtxPack::create` has improved, and we
///   can provide a detailed description as to why construction failed (without
///   writing much extra code), I'm starting to think that we should replace
///   this with a simple function that either returns a fully-constructed
///   GrackleCtxPack or aborts the program
///   - notably, this choice means that this gets rid of test-skipping
///     functionallity.
///   - I think the fact that this solution would abort the program rather than
///     gracefully reporting a problem is probably ok... Since error-reporting
///     is more robust, it would probably be straight-forward to write
///     alternative functions with specialized handling on a case-by-case basis
///     (e.g. for tests that check Grackle's behavior for deliberately invalid
///     configuration options)
/// - I'm inclined to stick with this macro for now, until the better handling
///   approach becomes a little more obvious...
#define GRTest_MAKE_CTX_PACK(lhs, conf)                                        \
  /* GRTest_CONCAT_TOK(rslt, __LINE__) is the temporary variable's name */     \
  std::pair<::grtest::GrackleCtxPack, ::grtest::Status> GRTest_CONCAT_TOKEN_(  \
      rslt, __LINE__) = ::grtest::GrackleCtxPack::create(conf);                \
                                                                               \
  /* Check whether construction succeeded */                                   \
  if (GRTest_CONCAT_TOKEN_(rslt, __LINE__).second.is_missing_std_file()) {     \
    GTEST_SKIP() << GRTest_CONCAT_TOKEN_(rslt, __LINE__).second.to_string();   \
  } else if (GRTest_CONCAT_TOKEN_(rslt, __LINE__).second.is_err()) {           \
    FAIL() << GRTest_CONCAT_TOKEN_(rslt, __LINE__).second.to_string()          \
           << "    occurred for:\n"                                            \
           << conf.stringify(false, "      ");                                 \
  }                                                                            \
                                                                               \
  /* store the result in lhs */                                                \
  lhs = std::move(GRTest_CONCAT_TOKEN_(rslt, __LINE__).first)

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
    GRTest_MAKE_CTX_PACK(pack, conf);
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
    ParamConf conf = GetParam();
    GRTest_MAKE_CTX_PACK(pack, conf);
  }

  GrackleCtxPack pack;
};

}  // namespace grtest

#endif  // GRTEST_FIXTURE
