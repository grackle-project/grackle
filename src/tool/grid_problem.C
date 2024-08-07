#include <array>
#include <cstdio>  // std::fprintf, stderr
#include <cstdlib> // std::exit
#include <regex>
#include <utility> // std::pair

#include "FieldData.h"
#include "grid_problem.h"
#include "view.h"

//-------------------------------------------------------------------

/// helper function that tries to parse a sequence
///
/// @returns The const char* component is a nullptr if there aren't any issues.
///   Otherwise, it specifies an error message.
static std::pair<seq_variant, const char*> parse_seq_(std::string s) {

  // create a dummy value to return if there are any issues
  seq_variant dummy_seq{valseq::Single(-1.0)};

  const std::regex flt_regex("^" FLT_PATTERN "$");

  if (std::regex_match(s, flt_regex)) {
    double val = std::stod(s);
    seq_variant seq{valseq::Single(val)};
    return { seq, nullptr };

  } else if (starts_with(s, "geomspace")) { // "geomspace;lo;hi;num"
    const char* pattern = (
        "^([a-zA-Z]+);(" FLT_PATTERN ");(" FLT_PATTERN ");([-+]?\\d+)$");
    std::regex geomspace_regex(pattern);
    std::smatch match;
    if (std::regex_match(s, match, geomspace_regex)) {
      // match 0 is full string
      // match 1 is "geomspace"

      // match 2 is the lo string
      double lo = std::stod(match[2].str());
      // matches 3-5 are ALWAYS various parts of float (even if we the string
      // doesn't have characters that match every subgroup)

      // match 6 is the hi string
      double hi = std::stod(match[6].str());
      // matches 7-9 are ALWAYS various parts of a float

      // match 10 is always the number of entries
      int num = std::stoi(match[10].str());
      //printf("geomspace;%g;%g;%d\n", lo,hi,num);
      seq_variant seq{valseq::Geomspace(lo,hi,num)};
      return {seq, nullptr};
    } else {
      return {dummy_seq,
              "gemospace must be formatted as geomspace;{lo};{hi};{num}"};
    }
  } else {
    return {dummy_seq,
            "sequence must be a lone value or a geomspace sequence"};
  }
}

static const std::regex grid_ax_regex("^-Sgrid.ax([0-9]+)=([a-zA-Z_]+);(.*)$");

static std::optional<scenario::QuantityKind> try_get_quantity_(std::string s) {
  if (s == "temperature")  return {scenario::QuantityKind::temperature};
  if (s == "mass_density") return {scenario::QuantityKind::mass_density};
  if (s == "metallicity")  return {scenario::QuantityKind::metallicity};
  return {};
}


bool scenario::try_parse_cli_grid_component(const char* arg,
                                            scenario::CliGridSpec& grid_spec)
{
  using namespace scenario;

  // if we ever add more scenario-kinds, we will use the following line
  //    if (!starts_with("-Sgrid")) return false;
  // in place of the next 2 if/elif cases
  if (!starts_with(arg, "-S")) {
    return false;
  } else if (!starts_with(arg, "-Sgrid")) {
    // if we ever add more scenario-kinds, 
    std::fprintf(
        stderr,
        "the argument %s specifies an argument belonging to an unknown "
        "kind of scenario. At this time, \"grid\" is the only known "
        "scenario kind",
        arg);
    std::exit(1);
  }

  // here we move onto actual parsing
  std::string_view ax_prefix("-Sgrid.ax");
  if (starts_with(arg, ax_prefix)) {
    std::string arg_str(arg);
    std::smatch match;
    if (std::regex_match(arg_str, match, grid_ax_regex)) {
      // match 0 is the full match

      // match 1 is the axis digit
      int axis_index = std::stoi(match[1].str());
      if ((axis_index < 0) || (axis_index > 2)) {
        std::fprintf(stderr,
                     "Error parsing `%s`. It must provide ax0, ax1, or ax2.\n",
                     arg);
        std::exit(1);
      } else if (grid_spec.ax_props[axis_index].has_value()) {
        std::fprintf(stderr,
                     "`%s` is the second argument to specify -Sgrid.ax%d\n",
                     arg, axis_index);
        std::exit(1);
      }

      // match 2 is the quantity-name
      std::optional<QuantityKind> quantity = try_get_quantity_(match[2].str());
      if (!quantity.has_value()) {
        std::fprintf(stderr,
                     "Error parsing `%s`. Invalid quantity name.\n", arg);
        std::exit(1);
      }

      // match 3 is the full sequence string
      std::string seq_str = match[3].str();
      std::pair<seq_variant, const char*> seq_rslt = parse_seq_(seq_str);
      if (seq_rslt.second != nullptr) {
        std::fprintf(stderr,
                     "Error parsing the sequence part of `%s`, or `%s`. %s\n",
                     arg, seq_str.c_str(), seq_rslt.second);
        std::exit(1);
      }

      AxisSpec ax_spec{quantity.value(), seq_rslt.first};
      grid_spec.ax_props[axis_index] = std::optional<AxisSpec>{ax_spec};

    } else {
      fprintf(
          stderr,
          "The argument `%s` is invalid. It looks like it should have the "
          "form `-Sgrid.ax<I>=<quantity>;<seq>`\n",
          arg);
      std::exit(1);
    }
  }

  // if we make it this far, then we have successfully parsed everything
  grid_spec.any_specified = true;
  return true;
}


//-------------------------------------------------------------------
// from here on out, we implement logic for actually setting up a grid

// the following should come from inside grackle
#define MU_METAL 16

static std::array<int,3> get_shape_(const scenario::CliGridSpec& scenario) {
  std::array<int, 3> shape;
  for (int i = 0; i < 3; i++) {
    if (scenario.ax_props[i].has_value()) {
      shape[i] = std::visit([](auto&& seq)->int { return seq.num; },
                            scenario.ax_props[i].value().seq);
    } else {
      shape[i] = 1;
    }
  }
  return shape;
}


static void initialize_View_from_ax_props_(
    View<gr_float> view, scenario::QuantityKind quantity,
    const scenario::CliGridSpec& scenario, double SolarMetalFractionByMass)
{
  bool complete = false;
  for (int ax_ind = 0; ax_ind < 3; ax_ind++ ) {
    if (!scenario.ax_props[ax_ind].has_value()) {
      break;
    } else if (scenario.ax_props[ax_ind]->quantity != quantity){
      continue;
    } else if (complete) {
      GRCLI_ERROR("a quantity was specifed for more than 1 axis");
    }
    
    double factor = (quantity != scenario::QuantityKind::metallicity)
                    ? SolarMetalFractionByMass : 1.0; 


    // define a lambda function to use for initialization. The argument is
    // templated on the kind of sequence
    auto fn = [ax_ind, factor, view](auto&& seq)-> void {

      for (int k = 0; k < view.extent(2); k++) {
        for (int j = 0; j < view.extent(1); j++) {
          for (int i = 0; i < view.extent(0); i++) {
            int seq_index = ((ax_ind == 0) * i + (ax_ind == 1) * j +
                             (ax_ind == 2) * k);
            view(i,j,k) = gr_float(seq.get(seq_index) * factor);
          }
        }
      }

    };

    std::visit(fn, scenario.ax_props[ax_ind].value().seq);
    complete = true;
  }

  GRCLI_REQUIRE(complete, "Unable to finish initializing view");
}

static void initialize_grid_(View<gr_float> target_temperature,
                             View<gr_float> target_density,
                             View<gr_float> target_metal_density,
                             FieldData& fields,
                             const chemistry_data& my_chem,
                             const code_units& initial_units) {

  View<gr_float> density       = fields.view("density");
  View<gr_float> eint          = fields.view("internal_energy");
  View<gr_float> metal_density = fields.view("metal_density");

  // now we load in species fields
  View<gr_float> HI_density    = fields.view("HI_density");
  View<gr_float> HII_density   = fields.view("HII_density");
  View<gr_float> HeI_density   = fields.view("HeI_density");
  View<gr_float> HeII_density  = fields.view("HeII_density");
  View<gr_float> HeIII_density = fields.view("HeIII_density");
  View<gr_float> e_density     = fields.view("e_density");

  View<gr_float> HM_density   = fields.view("HM_density");
  View<gr_float> H2I_density  = fields.view("H2I_density");
  View<gr_float> H2II_density = fields.view("H2II_density");

  View<gr_float> DI_density  = fields.view("DI_density");
  View<gr_float> DII_density = fields.view("DII_density");
  View<gr_float> HDI_density = fields.view("HDI_density");

  // we will use the nominal gamma everywhere (in the future, we can revisit
  // this and and take a more self-consistent gamma when you have molecular
  // Hydrogen)
  double gm1 = (my_chem.Gamma - 1.0);

  code_units tmp = initial_units; // make a copy to avoid a cast
  double eint_to_Tdivmu_factor = gm1 * get_temperature_units(&tmp);
  double Tdivmu_to_eint_factor = 1.0 / eint_to_Tdivmu_factor;

  const double tiny_number = 1e-10;
  const std::array<int, 3> shape = fields.grid_dimensions();
  // this borrows a lot from the Enzo-E test-problem

  for (int k = 0; k < shape[2]; k++) {
    for (int j = 0; j < shape[1]; j++) {
      for (int i = 0; i < shape[0]; i++) {

        double cur_density = target_density(i,j,k);
        double cur_metal_density = target_metal_density(i,j,k);

        density(i,j,k) = cur_density;
        metal_density(i,j,k) = cur_metal_density;

        if (my_chem.primordial_chemistry > 0){
          HI_density(i,j,k)    = cur_density * my_chem.HydrogenFractionByMass;
          HII_density(i,j,k)   = cur_density * tiny_number;
          // not going to worry about it right now, but it may make sense to
          // ignore contributions from metallicity
          HeI_density(i,j,k)   = cur_density * (1.0 - my_chem.HydrogenFractionByMass);
          HeII_density(i,j,k)  = cur_density * tiny_number;
          HeIII_density(i,j,k) = cur_density * tiny_number;
          e_density(i,j,k)     = cur_density * tiny_number;
        }
        if (my_chem.primordial_chemistry > 1){
          HM_density(i,j,k)    = cur_density * tiny_number;
          H2I_density(i,j,k)   = cur_density * tiny_number;
          H2II_density(i,j,k)  = cur_density * tiny_number;
        }
        if (my_chem.primordial_chemistry > 2){
          DI_density(i,j,k)    = cur_density * my_chem.DeuteriumToHydrogenRatio;
          DII_density(i,j,k)   = cur_density * tiny_number;
          HDI_density(i,j,k)   = cur_density * tiny_number;
        }

        // this is a natural place for improvements in the future:
        // -> 
        gr_float mu;
        if (my_chem.primordial_chemistry == 0) {
          mu = 0.6;
        } else {
          // we have already effectively assumed that everything is
          // neutral-atomic

          gr_float tmp = 
            ( (e_density(i,j,k) + HI_density(i,j,k) + HII_density(i,j,k)) +
              0.25 * (HeI_density(i,j,k) + HeII_density(i,j,k) +
                      HeIII_density(i,j,k) ) );

          if (my_chem.primordial_chemistry > 1) {
            tmp += HM_density(i,j,k) + 0.5 * (H2I_density(i,j,k) +
                                              H2II_density(i,j,k));
          }
          if (my_chem.primordial_chemistry > 2) {
            tmp += (0.5 * (DI_density(i,j,k) + DII_density(i,j,k))
                    + HDI_density(i,j,k)/3.0);
          }

          if (my_chem.metal_cooling == 1) tmp += cur_metal_density / MU_METAL;

          mu = density(i,j,k) / tmp;
        }

        // solve for eint from
        eint(i,j,k) = gr_float(
            Tdivmu_to_eint_factor * target_temperature(i,j,k) / mu);
      }
    }
  }

}

FieldData scenario::initialize_grid(const chemistry_data& my_chem,
                                    const code_units& initial_units,
                                    const scenario::CliGridSpec& scenario)
{

  // we will simply assume that my_chem is already initialized!

  // some basic checks!
  if (!scenario.ax_props[0].has_value()){
    GRCLI_ERROR("properties must always be provided for ax0");
  } else if (scenario.ax_props[2].has_value() && 
             !scenario.ax_props[1].has_value()) {
    GRCLI_ERROR("properties can't be provided for ax2 when ax1 isn't provided");
  }

  // get the shape of the fields
  std::array<int,3> shape = get_shape_(scenario);
  std::size_t size = shape[0] * shape[1] * shape[2];

  // create the FieldData
  FieldData fields(shape, my_chem);

  // we will now allocate arrays for each of the "axes" values
  // -> we could get a little more clever about this (we don't need a full 3D
  //    array for temperature)
  std::vector<gr_float> temperature_vec(size);
  View<gr_float> temperature(temperature_vec.data(), shape,
                             ContiguousLayout::Fortran);
  initialize_View_from_ax_props_(temperature,
                                 scenario::QuantityKind::temperature,
                                 scenario, my_chem.SolarMetalFractionByMass);

  std::vector<gr_float> density_vec(size);
  View<gr_float> density(density_vec.data(), shape, ContiguousLayout::Fortran); 
  initialize_View_from_ax_props_(density,
                                 scenario::QuantityKind::mass_density,
                                 scenario, my_chem.SolarMetalFractionByMass);

  std::vector<gr_float> metal_density_vec(size);
  View<gr_float> metal_density(metal_density_vec.data(), shape,
                               ContiguousLayout::Fortran); 
  initialize_View_from_ax_props_(metal_density,
                                 scenario::QuantityKind::metallicity,
                                 scenario, my_chem.SolarMetalFractionByMass);


  GRCLI_REQUIRE(my_chem.metal_cooling == 1,
      "we currently require metal cooling (can relax this in the future");

  initialize_grid_(temperature, density, metal_density, fields,
                   my_chem, initial_units);
  return fields;
}
