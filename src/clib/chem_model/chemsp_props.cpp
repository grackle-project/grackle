//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// implements logic pertaining to chemical species properties
///
//===----------------------------------------------------------------------===//

#include "chem_model/chemsp_props.hpp"

namespace GRIMPL_NAMESPACE_DECL {

Expected<ChemSpeciesProps, Error> ChemSpeciesProps::from_formula(
    std::string_view formula, const NuclideModel& model) {
  Expected<parse::FormulaComps, parse::Err> parse_rslt =
      parse::formula_comps(formula);
  if (!parse_rslt.has_value()) {
    return Unexpected(Error::parse(formula, parse_rslt.error())
                          .context_literal("problem parsing formula"));
  }
  return ChemSpeciesProps::from_formula(parse_rslt.value(), model);
}

Expected<ChemSpeciesProps, Error> ChemSpeciesProps::from_formula(
    const parse::FormulaComps& comps, const NuclideModel& model) {
  ChemSpeciesProps out;
  out.mass_factor = 0.0;
  int n_protons = 0;
  for (std::size_t i = 0; i < comps.symbol_subscript_pairs.size(); i++) {
    auto [name, coef] = comps.symbol_subscript_pairs[i];

    // lookup the id associated with the nuclide specified by name
    std::optional<uint16_t> maybe_id = model.symbol_map().find(name);
    if (!maybe_id.has_value()) {
      std::string name_str(name);
      return Unexpected(Error::msgf("\"%s\" is not a known chemical species",
                                    name_str.c_str()));
    }
    uint16_t nuclide_id = maybe_id.value();

    // retrieve the nuclide's properties
    const NuclideProp& nuclide_prop = model.get(nuclide_id);

    // now, let's update total number of protons and update out
    n_protons += static_cast<int>(coef) * nuclide_prop.n_proton;
    out.mass_factor += static_cast<double>(coef) * nuclide_prop.mass_factor;
    if (!out.nuclide_id_coef_pairs.try_emplace_back(nuclide_id, coef)) {
      // shouldn't be reachable since that max capacity should be identical
      return Unexpected(Error::msgf(
          "nuclide_id_coef_pairs doesn't have adequate capacity to hold "
          "data for nuclide symbol %zu",
          i));
    }
  }
  out.elementary_charge = comps.roman_numeral_val - 1;

  // perform one last sanity check (perhaps we should move it to a unit test)
  if ((out.elementary_charge - n_protons) > 0) {
    return Unexpected(Error::msgf(
        "species can't have a charge of +%d (it only contains %d protons)",
        out.elementary_charge, n_protons));
  }
  return out;
}

}  // namespace GRIMPL_NAMESPACE_DECL