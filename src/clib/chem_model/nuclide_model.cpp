//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Defines NuclideModel
///
//===----------------------------------------------------------------------===//

#include <utility>  // std::pair

#include "nuclide_model.hpp"
#include "../phys_constants.hpp"
#include "../support/FrozenKeyIdxBiMap.hpp"
#include "../support/status_reporting.hpp"

namespace GRIMPL_NAMESPACE_DECL {

namespace nuclide_detail {

using Pair = std::pair<const char*, NuclideProp>;
constexpr Pair pairs[] = {
    // todo: we may want to consider changing the recorded mass of "H"
    // -> we probably want to use the mass for protium (i.e. the Hydrogen
    //    isotope without neutrons)
    // -> this may mess with python answer tests
    {"H", NuclideProp(1, true, true, 1.00794, mass_factor::H)},
    {"D", NuclideProp(1, true, true, 2.0141017778, mass_factor::D)},
    {"He", NuclideProp(2, true, false, 4.002602, mass_factor::He)},
    {"Li", NuclideProp(3, false, false, 6.941, MostCommonMassNumber{7})},
    {"Be", NuclideProp(4, false, false, 9.012182, MostCommonMassNumber{9})},
    {"B", NuclideProp(5, false, false, 10.811, MostCommonMassNumber{11})},
    {"C", NuclideProp(6, false, false, 12.0107, mass_factor::C)},
    {"N", NuclideProp(7, false, false, 14.0067, MostCommonMassNumber{14})},
    {"O", NuclideProp(8, false, false, 15.9994, mass_factor::O)},
    {"F", NuclideProp(9, false, false, 18.9984032, MostCommonMassNumber{19})},
    {"Ne", NuclideProp(10, false, false, 20.1797, MostCommonMassNumber{20})},
    {"Na", NuclideProp(11, false, false, 22.98977, MostCommonMassNumber{23})},
    {"Mg", NuclideProp(12, false, false, 24.305, mass_factor::Mg)},
    {"Al", NuclideProp(13, false, false, 26.981538, mass_factor::Al)},
    {"Si", NuclideProp(14, false, false, 28.0855, mass_factor::Si)},
    {"P", NuclideProp(15, false, false, 30.973761, MostCommonMassNumber{31})},
    {"S", NuclideProp(16, false, false, 32.065, mass_factor::S)},
    {"Cl", NuclideProp(17, false, false, 35.453, MostCommonMassNumber{35})},
    {"Ar", NuclideProp(18, false, false, 39.948, MostCommonMassNumber{40})},
    {"K", NuclideProp(19, false, false, 39.0983, MostCommonMassNumber{39})},
    {"Ca", NuclideProp(20, false, false, 40.078, MostCommonMassNumber{40})},
    {"Sc", NuclideProp(21, false, false, 44.95591, MostCommonMassNumber{45})},
    {"Ti", NuclideProp(22, false, false, 47.867, MostCommonMassNumber{48})},
    {"V", NuclideProp(23, false, false, 50.9415, MostCommonMassNumber{51})},
    {"Cr", NuclideProp(24, false, false, 51.9961, MostCommonMassNumber{52})},
    {"Mn", NuclideProp(25, false, false, 54.938049, MostCommonMassNumber{55})},
    {"Fe", NuclideProp(26, false, false, 55.845, mass_factor::Fe)},
    {"Co", NuclideProp(27, false, false, 58.9332, MostCommonMassNumber{59})},
    {"Ni", NuclideProp(28, false, false, 58.6934, MostCommonMassNumber{58})},
    {"Cu", NuclideProp(29, false, false, 63.546, MostCommonMassNumber{63})},
    {"Zn", NuclideProp(30, false, false, 65.409, MostCommonMassNumber{64})},
};

static constexpr int N_ENTRIES = static_cast<int>(sizeof(pairs) / sizeof(Pair));

}  // namespace nuclide_detail

NuclideModel::NuclideModel() {
  // initialize symbol_map_
  const char* symbols[nuclide_detail::N_ENTRIES];
  for (int i = 0; i < nuclide_detail::N_ENTRIES; i++) {
    symbols[i] = nuclide_detail::pairs[i].first;
  }
  symbol_map_ = FrozenKeyIdxBiMap::create(symbols, nuclide_detail::N_ENTRIES,
                                          BiMapMode::REFS_KEYDATA);
  if (!symbol_map_.is_ok()) {  // this is an internal programming error!
    GR_INTERNAL_ERROR("issue initializing symbol_map_");
  }

  // initialize props_
  props_.reserve(nuclide_detail::N_ENTRIES);
  for (int i = 0; i < nuclide_detail::N_ENTRIES; i++) {
    props_.push_back(nuclide_detail::pairs[i].second);
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL
