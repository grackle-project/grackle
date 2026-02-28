//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Test the API for dynamically querying `chemistry_data`
///
//===----------------------------------------------------------------------===//

#include "grackle.h"
#include "grtestutils/iterator_adaptor.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>

TEST(DynamicChemParam, EmptyStrNumParams) {
  EXPECT_EQ(grackle_num_params(""), 0);
}

TEST(DynamicChemParam, NullptrNumParams) {
  EXPECT_EQ(grackle_num_params(nullptr), 0);
}

TEST(DynamicChemParam, NotATypeNumParams) {
  EXPECT_EQ(grackle_num_params("NotARealType"), 0);
}

TEST(DynamicChemParam, ParamNameIndexZero) {
  EXPECT_NE(param_name_int(0u), nullptr);
  EXPECT_NE(param_name_double(0u), nullptr);
  EXPECT_NE(param_name_string(0u), nullptr);
}

TEST(DynamicChemParam, ParamNameIndexMaxNumber) {
  EXPECT_NE(param_name_int(grackle_num_params("int") - 1u), nullptr);
  EXPECT_NE(param_name_double(grackle_num_params("double") - 1u), nullptr);
  EXPECT_NE(param_name_string(grackle_num_params("string") - 1u), nullptr);
}

TEST(DynamicChemParam, ParamNameIndexPassCount) {
  EXPECT_EQ(param_name_int(grackle_num_params("int")), nullptr);
  EXPECT_EQ(param_name_double(grackle_num_params("double")), nullptr);
  EXPECT_EQ(param_name_string(grackle_num_params("string")), nullptr);
}

TEST(DynamicChemParam, ParamNameIndexPassCountPlus1) {
  EXPECT_EQ(param_name_int(grackle_num_params("int") + 1u), nullptr);
  EXPECT_EQ(param_name_double(grackle_num_params("double") + 1u), nullptr);
  EXPECT_EQ(param_name_string(grackle_num_params("string") + 1u), nullptr);
}

TEST(DynamicChemParam, ParamNameMaximum) {
  unsigned i = std::numeric_limits<unsigned>::max();

  EXPECT_EQ(param_name_int(i), nullptr);
  EXPECT_EQ(param_name_double(i), nullptr);
  EXPECT_EQ(param_name_string(i), nullptr);
}

static constexpr std::array<grtest::ChemParamType, 3> known_chem_types = {
    grtest::ChemParamType::INT, grtest::ChemParamType::DOUBLE,
    grtest::ChemParamType::STRING};

TEST(DynamicChemParam, UniqueNames) {
  // validate that all parameter names are unique
  std::set<std::string> name_set;

  for (const grtest::ChemParamType type : known_chem_types) {
    // iterate over all parameter names of the current type
    for (const auto [param_name, _] : grtest::ChemParamRange(type)) {
      ASSERT_TRUE(name_set.insert(param_name).second)
          << "the name, \"" << param_name << "\" appears more than once";
    }
  }
}

TEST(DynamicChemParam, AccessNullChemistryData) {
  EXPECT_EQ(local_chemistry_data_access_int(nullptr, "use_grackle"), nullptr);
  EXPECT_EQ(local_chemistry_data_access_double(nullptr, "Gamma"), nullptr);
  EXPECT_EQ(local_chemistry_data_access_string(nullptr, "grackle_data_file"),
            nullptr);
}

TEST(DynamicChemParam, AccessEmptyString) {
  chemistry_data my_chemistry;
  ASSERT_EQ(local_initialize_chemistry_parameters(&my_chemistry), GR_SUCCESS);

  const char* name = "";

  EXPECT_EQ(local_chemistry_data_access_int(&my_chemistry, name), nullptr);
  EXPECT_EQ(local_chemistry_data_access_double(&my_chemistry, name), nullptr);
  EXPECT_EQ(local_chemistry_data_access_string(&my_chemistry, name), nullptr);
}

TEST(DynamicChemParam, AccessNotAParam) {
  chemistry_data my_chemistry;
  ASSERT_EQ(local_initialize_chemistry_parameters(&my_chemistry), GR_SUCCESS);

  const char* name = "NotARealParameter";

  EXPECT_EQ(local_chemistry_data_access_int(&my_chemistry, name), nullptr);
  EXPECT_EQ(local_chemistry_data_access_double(&my_chemistry, name), nullptr);
  EXPECT_EQ(local_chemistry_data_access_string(&my_chemistry, name), nullptr);
}

static void* access_ptr_(chemistry_data* my_chem, grtest::ChemParamType type,
                         const char* param_name) {
  switch (type) {
    case grtest::ChemParamType::INT: {
      return (void*)local_chemistry_data_access_int(my_chem, param_name);
    }
    case grtest::ChemParamType::DOUBLE: {
      return (void*)local_chemistry_data_access_double(my_chem, param_name);
    }
    case grtest::ChemParamType::STRING: {
      return (void*)local_chemistry_data_access_string(my_chem, param_name);
    }
    default:
      return nullptr;
  }
}

TEST(DynamicChemParam, UniqueAddresses) {
  chemistry_data my_chem;
  ASSERT_EQ(local_initialize_chemistry_parameters(&my_chem), GR_SUCCESS);

  // we will validate that all parameter names are successfully access addresses
  // with different names
  const std::map<grtest::ChemParamType, std::string> map({
      {grtest::ChemParamType::INT, "int"},
      {grtest::ChemParamType::DOUBLE, "double"},
      {grtest::ChemParamType::STRING, "string"},
  });
  std::set<uintptr_t> address_set;

  for (const grtest::ChemParamType type : known_chem_types) {
    // iterate over all parameter names of the current type
    for (const auto [param_name, _] : grtest::ChemParamRange(type)) {
      void* ptr = access_ptr_(&my_chem, type, param_name.c_str());
      EXPECT_NE(ptr, nullptr)
          << "local_chemistry_data_access_" << map.at(type)
          << " returned nullptr for the " << map.at(type) << "-typed "
          << "parameter, \"" << param_name << '"';

      uintptr_t address = reinterpret_cast<uintptr_t>(ptr);
      ASSERT_TRUE(address_set.insert(address).second)
          << "address associated with, \"" << param_name << "\" appears more "
          << "than once";
    }
  }
}

TEST(DynamicChemParam, AccessWrongType) {
  chemistry_data my_chem;
  ASSERT_EQ(local_initialize_chemistry_parameters(&my_chem), GR_SUCCESS);

  for (const grtest::ChemParamType type : known_chem_types) {
    // iterate over all parameter names of the current type
    for (const auto [param_name, _] : grtest::ChemParamRange(type)) {
      if (type != grtest::ChemParamType::INT) {
        int* ptr =
            local_chemistry_data_access_int(&my_chem, param_name.c_str());
        EXPECT_EQ(ptr, nullptr)
            << "local_chemistry_data_access_int should "
            << "return a nullptr when called to access a non-integer "
            << "parameter, like \"" << param_name << '"';
      }

      if (type != grtest::ChemParamType::DOUBLE) {
        double* ptr =
            local_chemistry_data_access_double(&my_chem, param_name.c_str());
        EXPECT_EQ(ptr, nullptr)
            << "local_chemistry_data_access_double should "
            << "return a nullptr when called to access a non-double "
            << "parameter, like \"" << param_name << '"';
      }

      if (type != grtest::ChemParamType::STRING) {
        char** ptr =
            local_chemistry_data_access_string(&my_chem, param_name.c_str());
        EXPECT_EQ(ptr, nullptr)
            << "local_chemistry_data_access_string should "
            << "return a nullptr when called to access a non-string "
            << "parameter, like \"" << param_name << '"';
      }
    }
  }
}
