//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Check correctness of the index_helper logic
///
//===----------------------------------------------------------------------===//

#include <optional>
#include <ostream>
#include <vector>

#include <gtest/gtest.h>

#include "grackle.h"
#include "support/index_helper.hpp"
#include "support/status_reporting.hpp"

#include "grtestutils/view.hpp"
#include "grtestutils/googletest/check_allclose.hpp"

namespace {  // stuff inside an anonymous namespace is local to this file

enum class IndexingStrat { FLAT, MULTIDIM };

const char* to_str_(const IndexingStrat& strat) {
  switch (strat) {
    case IndexingStrat::FLAT:
      return "IndexingStrat::FLAT";
    case IndexingStrat::MULTIDIM:
      return "IndexingStrat::MULTIDIM";
    default:
      GR_INTERNAL_UNREACHABLE_ERROR();
  }
}

// teach GoogleTest how to print IndexingStrat
void PrintTo(const IndexingStrat& s, std::ostream* os) { *os << to_str_(s); }

/// @brief A class that represents an index region
///
/// This only exists to help write tests in a relatively concise manner
class IndexRegion {
  friend std::vector<double> fill_mask(const IndexRegion& idx_region,
                                       IndexingStrat s);

  // the attributes
  grtest::IdxMapping<grtest::DataLayout::LEFT> idx_mapping_;
  std::vector<int> start_;
  std::vector<int> stop_;

public:
  IndexRegion(grtest::IdxMapping<grtest::DataLayout::LEFT> idx_mapping,
              std::vector<int> start, std::vector<int> stop)
      : idx_mapping_(idx_mapping), start_(start), stop_(stop) {
    int rank = idx_mapping_.rank();
    GR_INTERNAL_REQUIRE(rank >= 1, "non-positive rank is forbidden");
    GR_INTERNAL_REQUIRE(start_.size() == static_cast<std::size_t>(rank),
                        "the number of args in start arg must match rank");
    GR_INTERNAL_REQUIRE(stop_.size() == static_cast<std::size_t>(rank),
                        "the number of args in stop arg must match rank");
  }
};

/// return a mask representative of the index region
///
/// The mask has an entry for each location in the idx_region. Each entry is
/// initialized to 0 and incremented by 1 when the loop visits the location
std::vector<double> fill_mask(const IndexRegion& idx_region, IndexingStrat s) {
  int rank = idx_region.idx_mapping_.rank();

  // set up index_helper:
  GRIMPL_NS::IndexHelper index_helper;
  {
    std::vector<int> start = idx_region.start_;
    std::vector<int> end(rank);
    std::vector<int> extent(rank);
    for (int i = 0; i < rank; i++) {
      end[i] = idx_region.stop_[i] - 1;
      extent[i] = idx_region.idx_mapping_.extents()[i];
    }

    grackle_field_data tmp;
    gr_initialize_field_data(&tmp);
    tmp.grid_rank = rank;
    tmp.grid_dimension = extent.data();
    tmp.grid_start = start.data();
    tmp.grid_end = end.data();

    index_helper = GRIMPL_NS::build_index_helper_(&tmp);
  }
  int outer_idx_size = index_helper.outer_ind_size;

  // construct the mask and then fill it
  std::vector<double> out(idx_region.idx_mapping_.n_elements());
  switch (s) {
    case IndexingStrat::FLAT: {
      for (int outer_idx = 0; outer_idx < outer_idx_size; outer_idx++) {
        GRIMPL_NS::FieldFlatIndexRange flat_idx_range =
            inner_flat_range_(outer_idx, &index_helper);
        for (int i = flat_idx_range.start; i <= flat_idx_range.end; i++) {
          out[i] += 1;
        }
      }
      break;
    }
    case IndexingStrat::MULTIDIM: {
      const grtest::IdxMapping<grtest::DataLayout::LEFT>& idx_mapping =
          idx_region.idx_mapping_;
      if (rank == 1) {
        for (int outer_idx = 0; outer_idx < outer_idx_size; outer_idx++) {
          GRIMPL_NS::IndexRange idx_range =
              GRIMPL_NS::make_idx_range_(outer_idx, &index_helper);
          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            out[i] += 1;
          }
        }
      } else if (rank == 2) {
        for (int outer_idx = 0; outer_idx < outer_idx_size; outer_idx++) {
          GRIMPL_NS::IndexRange idx_range =
              GRIMPL_NS::make_idx_range_(outer_idx, &index_helper);
          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            int flat_idx = idx_mapping(i, idx_range.j);
            out[flat_idx] += 1;
          }
        }
      } else if (rank == 3) {
        for (int outer_idx = 0; outer_idx < outer_idx_size; outer_idx++) {
          GRIMPL_NS::IndexRange idx_range =
              GRIMPL_NS::make_idx_range_(outer_idx, &index_helper);
          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            int flat_idx = idx_mapping(i, idx_range.j, idx_range.k);
            out[flat_idx] += 1;
          }
        }
      } else {
        GR_INTERNAL_UNREACHABLE_ERROR();
      }
      break;
    }
    default:
      GR_INTERNAL_UNREACHABLE_ERROR();
  }

  return out;
}

testing::AssertionResult run_full_test(std::vector<int> shape,
                                       std::vector<int> start,
                                       std::vector<int> stop,
                                       IndexingStrat strat,
                                       const double* expected) {
  std::optional<grtest::IdxMapping<grtest::DataLayout::LEFT>> maybe_idx_map =
      grtest::IdxMapping<grtest::DataLayout::LEFT>::try_create(shape.size(),
                                                               shape.data());
  if (!maybe_idx_map.has_value()) {
    GR_INTERNAL_ERROR("bad shape arg");
  }
  grtest::IdxMapping<grtest::DataLayout::LEFT> idx_map = maybe_idx_map.value();
  IndexRegion region(idx_map, start, stop);
  std::vector<double> access_locs = fill_mask(region, strat);

  return check_array_equal(access_locs.data(), expected, idx_map);
}

}  // anonymous namespace

class IndexHelperTest : public testing::TestWithParam<IndexingStrat> {};

TEST_P(IndexHelperTest, Simple1D) {
  std::vector<int> shape = {5};
  std::vector<int> start = {0};
  std::vector<int> stop = shape;
  const double expected[] = {1, 1, 1, 1, 1};

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipLead1D) {
  std::vector<int> shape = {5};
  std::vector<int> start = {3};
  std::vector<int> stop = shape;
  const double expected[] = {0, 0, 0, 1, 1};

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipTrail1D) {
  std::vector<int> shape = {5};
  std::vector<int> start = {0};
  std::vector<int> stop = {4};
  const double expected[] = {1, 1, 1, 1, 0};

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipAllButOne1D) {
  std::vector<int> shape = {5};
  std::vector<int> start = {2};
  std::vector<int> stop = {3};
  const double expected[] = {0, 0, 1, 0, 0};

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, Simple2D) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {0, 0};
  std::vector<int> stop = shape;
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipLead2DFastAx) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {2, 0};
  std::vector<int> stop = shape;
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    0, 0, 1, 1, 1,
    0, 0, 1, 1, 1,
    0, 0, 1, 1, 1,
    0, 0, 1, 1, 1
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipLead2DSlowAx) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {0, 1};
  std::vector<int> stop = shape;
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    0, 0, 0, 0, 0,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipLead2DBothAx) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {1, 3};
  std::vector<int> stop = shape;
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 1, 1, 1, 1,
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipTrail2DFastAx) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {0, 0};
  std::vector<int> stop = {3, 4};
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    1, 1, 1, 0, 0,
    1, 1, 1, 0, 0,
    1, 1, 1, 0, 0,
    1, 1, 1, 0, 0
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipTrail2DSlowAx) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {0, 0};
  std::vector<int> stop = {5, 2};
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, SkipTrail2DBothAx) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {0, 0};
  std::vector<int> stop = {2, 3};
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 0,
    0, 0, 0, 0, 0
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, Complex2D) {
  std::vector<int> shape = {5, 4};
  std::vector<int> start = {1, 2};
  std::vector<int> stop = {3, 3};
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 1, 1, 0, 0,
    0, 0, 0, 0, 0
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, Simple3D) {
  std::vector<int> shape = {5, 4, 3};
  std::vector<int> start = {0, 0, 0};
  std::vector<int> stop = shape;
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,

    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,

    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

TEST_P(IndexHelperTest, Complex3D) {
  std::vector<int> shape = {5, 4, 3};
  std::vector<int> start = {1, 2, 1};
  std::vector<int> stop = {3, 3, 2};
  // clang-format off:  formatter knows nothing about shape
  const double expected[] = {
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,

    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 1, 1, 0, 0,
    0, 0, 0, 0, 0,

    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
  };
  // clang-format on

  EXPECT_TRUE(run_full_test(shape, start, stop, GetParam(), expected));
}

constexpr static const IndexingStrat strats_[] = {IndexingStrat::FLAT,
                                                  IndexingStrat::MULTIDIM};
INSTANTIATE_TEST_SUITE_P(, IndexHelperTest, testing::ValuesIn(strats_));