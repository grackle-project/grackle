#include <array>

#include <view.h>
#include "grid_problem.h"

// the following should come from inside grackle
#define MU_METAL 16

static std::array<int,3> get_shape_(const GridScenarioSpec& scenario) {
  std::array<int, 3> shape;
  for (int i = 0; i < 3; i++) {
    if (scenario.ax_props[i].has_value()) {
      shape[i] = std::visit([](auto&& seq)->int { return seq.num; },
                            scenario.ax_props[i].has_value());
    } else {
      shape[i] = 1;
    }
  }
  return shape;
}


static void initialize_View_from_ax_props_(View<gr_float> view,
                                           QuantityKind quantity,
                                           const GridScenarioSpec& scenario,
                                           double SolarMetalFractionByMass)
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
    
    double factor = (quantity != QuantityKind::metallicity)
                    ? SolarMetalFractionByMass : 1.0; 


    // define a lambda function to use for initialization. The argument is
    // templated on the kind of sequence
    auto fn = [ax_ind, factor, view](auto&& seq)-> void {

      for (int k = 0; k < view.shape[2]; k++) {
        for (int j = 0; j < view.shape[1]; j++) {
          for (int i = 0; i < view.shape[0]; i++) {
            int seq_index = ((ax_ind == 0) * i + (ax_ind == 1) * j +
                             (ax_ind == 2) * k);
            view(i,j,k) = gr_float(seq.get(seq_index) * factor);
          }
        }
      }

    };

    std::visit(fn, scenario.ax_props[i]->seq);
    complete = true;
  }

  GRCLI_REQUIRE(complete, "Unable to finish initializing view");
}


void initialize_grid_(View<gr_float> target_temperature,
                      View<gr_float> target_density,
                      View<gr_float> target_metal_density,
                      FieldData& fields,
                      const chemistry_data& my_chem,
                      const code_units& initial_units) {

  View<gr_float> density       = field.view("density");
  View<gr_float> eint          = field.view("internal_energy");
  View<gr_float> metal_density = field.view("metal_density");

  // now we load in species fields
  View<gr_float> HI_density    = field.view("HI_density");
  View<gr_float> HII_density   = field.view("HII_density");
  View<gr_float> HeI_density   = field.view("HeI_density");
  View<gr_float> HeII_density  = field.view("HeII_density");
  View<gr_float> HeIII_density = field.view("HeIII_density");
  View<gr_float> e_density     = field.view("e_density");

  View<gr_float> HM_density   = field.view("HM_density");
  View<gr_float> H2I_density  = field.view("H2I_density");
  View<gr_float> H2II_density = field.view("H2II_density");

  View<gr_float> DI_density  = field.view("DI_density");
  View<gr_float> DII_density = field.view("DII_density");
  View<gr_float> HDI_density = field.view("HDI_density");


  // get the conversion factor (make a copy my_units so we don't need to cast)
  code_units tmp = initial_units;
  double eint_to_Tdivmu_factor 

  const double tiny_number = 1e-10;
  const std::array<int, 3> shape = field.grid_dimensions();
  // this borrows a lot from the Enzo-E test-problem

  for (int k = 0; k < shape[2]; k++) {
    for (int j = 0; j < shape[1]; j++) {
      for (int i = 0; i < shape[0]; i++) {

        double cur_density = nominal_density(i,j,k);
        double cur_metal_density = target_metal_density(i,j,k);

        density(i,j,k) = cur_density;
        metal_density(i,j,k) = cur_metal_density;

        if (my_chem.primordial_chemistry > 0){
          HI_density(i,j,k)    = cur_density * my_chem.HydrogenFractionByMass;
          HII_density(i,j,k)   = cur_density * tiny_number;
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


        // we may want to support customization in the future...

        gr_float mu;
        if (my_chem.primordial_chemistry == 0) {
          mu = 0.6;
        } else {
          // we have already effectively assumed that everything is
          // neutral-atomic

          gr_float tmp = 
            ( (e_density(iz,iy,ix) + HI_density(iz,iy,ix) +
               HII_density(iz,iy,ix)) +
              0.25 * (HeI_density(iz,iy,ix) + HeII_density(iz,iy,ix) +
                      HeIII_density(iz,iy,ix) ) );

          if (my_chem.primordial_chemistry > 1) {
            tmp += HM_density(iz,iy,ix) + 0.5 * (H2I_density(iz,iy,ix) +
                                                 H2II_density(iz,iy,ix));
          }
          if (my_chem.primordial_chemistry > 2) {
            tmp += (0.5 * (DI_density(iz,iy,ix) + DII_density(iz,iy,ix))
                    + HDI_density(iz,iy,ix)/3.0);
          }

          if (my_chem.metal_cooling == 1) {
            tmp += cur_metal_density / MU_METAL;
          }

          mu = density(iz,iy,ix) / tmp;
        }

        // we will use the nominal gamma
        gr_float gm1 = (my_chem.Gamma - 1.0);

        get_temperature_units(my_units)

      }
    }
  }

}

grackle_field_data* initialize_grid(const chemistry_data& my_chem,
                                    const code_units& initial_units,
                                    const GridScenarioSpec& scenario) {
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

  grackle_field_data* fields = impl::allocate_and_init_gr_field_data
    (my_chem, rank, shape.data());

  // we will now allocate arrays for each of the "axes" values
  // -> we could get a little more clever about this (we don't need a full 3D
  //    array for temperature)
  std::vector<gr_float> temperature_vec(size);
  View<gr_float> temperature(temperature_vec.data(), shape,
                             ContiguousLayout::Fortran)
  initialize_View_from_ax_props_(temperature, QuantityKind::temperature,
                                 scenario, my_chem.SolarMetalFractionByMass);

  std::vector<gr_float> density_vec(size);
  View<gr_float> density(density_vec.data(), shape, ContiguousLayout::Fortran); 
  initialize_View_from_ax_props_(density, QuantityKind::mass_density,
                                 scenario, my_chem.SolarMetalFractionByMass);

  std::vector<gr_float> metal_density_vec(size);
  View<gr_float> metal_density(metal_density_vec.data(), shape,
                               ContiguousLayout::Fortran); 
  initialize_View_from_ax_props_(metal_density, QuantityKind::metal_density,
                                 scenario, my_chem.SolarMetalFractionByMass);


  GRCLI_REQUIRE(my_chem.metal_cooling == 1,
                "we currently require metal cooling (can relax this in the
                "future");

  gr_float* metal_density = new gr_float[size];




}

