//@% IF YOU ARE A NEW USER: you should go to the documentation page where we
//@% describe this example
//@%
//@% This is a minimal example of how to use Grackle:
//@% - we talk about sections of logic chunk-by-chunk in the documentation
//@% - lines starting with ``//@!`` are not rendered in the documentation.
//@% - lines of the form ``//@![[[ANCHOR:<description>]]]`` are generally used
//@%   as "anchors" to help us include snippets of text in the documentation
//@%
//@% Whenever you change the contents of this file, you should ensure that the
//@% the corresponding page of documentation still makes sense

#include <grackle.h>
#include <stdio.h>
#include <stdlib.h>

const char* path = "../../input/CloudyData_UVB=HM2012.h5";  // may need changing
typedef double Real;  // <- precision of our field data

int main(int argc, char** argv) {
  //@![[[ANCHOR:ScenarioSetup]]]
  const double densU = 1.67e-24, lenU = 3.086e21, timeU = 3.16e13;
  double energyU = (lenU*lenU)/(timeU*timeU);   // <- internal energy units
  int dims[3]  = {1, 1, 1};                     // <- shape of our fields
  Real rho[1]  = { (Real)(2.0e-24 / densU) };   // <- density field data
  Real eint[1] = { (Real)(2.0e12 / energyU) };  // <- internal energy field data

  // we want to apply radiative cooling over the following timestep
  double dt = 1.0e13 / timeU; // corresponds to ~0.32 Myr

  //@![[[ANCHOR:SanityChecks]]]
  if (sizeof(gr_float) != sizeof(Real))  abort();     // config sanity check!
  if (gr_check_consistency() != GR_SUCCESS)  abort(); // linking sanity check!

  //@![[[ANCHOR:SetupChemistrySolver]]]
  // set up the chemistry solver
  chemistry_data par; // holds runtime parameters
  if (local_initialize_chemistry_parameters(&par) != GR_SUCCESS)  abort();
  par.use_grackle  = 1;             par.with_radiative_cooling = 1;
  par.Gamma        = 5.0/3.0;       par.primordial_chemistry   = 0;
  par.UVbackground = 1;             par.grackle_data_file      = path;
  code_units u = {.comoving_coordinates=0, .a_units=1.0, .a_value=1.0,
                  .density_units=densU, .length_units=lenU, .time_units=timeU};
  chemistry_data_storage sto; // <- stores data tables & rate data
  if (local_initialize_chemistry_data(&par, &sto, &u) != GR_SUCCESS)  abort();

  //@![[[ANCHOR:SetupGrackleFieldData]]]
  grackle_field_data f; // <- used to communicate field data to solver
  if (gr_initialize_field_data(&f) != GR_SUCCESS)  abort();
  int start[3] = {0, 0, 0};
  int end[3] = {dims[0]-1, dims[1]-1, dims[2]-1};
  f.grid_dimension = dims; f.grid_start = start; f.grid_end = end;
  f.grid_rank = 3;         f.density    = rho;   f.internal_energy = eint;

  //@![[[ANCHOR:ApplySolver]]]
  // apply the solver
  if (local_solve_chemistry(&par, &sto, &u, &f, dt) != GR_SUCCESS) abort();
  printf("energy after timestep: %24.16e cm^2/s^2\n", eint[0]*energyU);

  //@![[[ANCHOR:Exit]]]
  return 0; // we skipped cleanup of solver (see local_free_chemistry_data)
}
