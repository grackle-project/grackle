// See LICENSE file for license and copyright information

/// @file solve_rate_cool_g-cpp.C
/// @brief Declares signature of solve_rate_cool_g

// This file was initially generated automatically during conversion of the
// solve_rate_cool_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_wrappers.hpp"
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "utils-cpp.hpp"

#include "solve_rate_cool_g-cpp.h"

/// Computes the timescale given by `ndens_Heq / (d ndens_Heq / d t)`
///
/// This is used to help compute the subcycle timestep when using a primordial
/// chemistry solver
///
/// @param my_chemistry holds a number of configuration parameters
/// @param my_rates holds assorted rate data. In this function, this is being
///    used to specify some the interpolation tables of some relevant reaction
///    rates (they are tabulated with respect to logT)
/// @param dlogtem Specifies the constant spacing shared by the relevant rate
///    interpolation tables
/// @param logTlininterp_buf Specifies the information related to the position
///    in the logT interpolations (for a number of chemistry zones)
/// @param k13, k22 1D arrays specifying the previously looked up, local values
///    of the k13 and k22 rates.
/// @param local_rho specifies the local (total) mass density
/// @param tgas, p2d, edot 1D arrays containing local values of temperature,
///    pressure divided by density, and the time derivative of internal energy.
/// @param i Specifies the index of the relevant zone in the 1D array. (**BE
///    AWARE:** this is a 0-based index)
///
/// @note
/// The `static` annotation indicates that this function is only visible to the
/// current translation unit
static double calc_Heq_div_dHeqdt_(
  const chemistry_data* my_chemistry,
  const chemistry_data_storage* my_rates,
  double dlogtem,
  const grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  const double* k13,
  const double* k22,
  double local_rho,
  const double* tgas,
  const double* p2d,
  const double* edot,
  int i
) {

  // Equilibrium value for H is:
  // Heq = (-1._DKIND / (4*k22)) * (k13 - sqrt(8 k13 k22 rho + k13^2))
  // We want to know dH_eq/dt.
  // - We can trivially get dH_eq/dT.
  // - We have de/dt.
  // - We need dT/de.
  //
  // T = (g-1)*p2d*utem/N; tgas == (g-1)(p2d*utem/N)
  // dH_eq / dt = (dH_eq/dT) * (dT/de) * (de/dt)
  // dH_eq / dT (see above; we can calculate the derivative here)
  // dT / de = utem * (gamma - 1._DKIND) / N == tgas / p2d
  // de / dt = edot
  // Now we use our estimate of dT/de to get the estimated
  // difference in the equilibrium
  double eqt2 = std::fmin(std::log(tgas[i]) + 0.1*dlogtem, logTlininterp_buf.t2[i]);
  double eqtdef = (eqt2 - logTlininterp_buf.t1[i])/(logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);
  double eqk222 = my_rates->k22[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k22[logTlininterp_buf.indixe[i]+1-1] -my_rates->k22[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double eqk132 = my_rates->k13[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k13[logTlininterp_buf.indixe[i]+1-1] -my_rates->k13[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double heq2 = (-1. / (4.*eqk222)) * (eqk132-
    std::sqrt(8.*eqk132*eqk222*
              my_chemistry->HydrogenFractionByMass*local_rho+
              std::pow(eqk132,2.)));

  double eqt1 = std::fmax(std::log(tgas[i]) - 0.1*dlogtem, logTlininterp_buf.t1[i]);
  eqtdef = (eqt1 - logTlininterp_buf.t1[i])/(logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);
  double eqk221 = my_rates->k22[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k22[logTlininterp_buf.indixe[i]+1-1] -my_rates->k22[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double eqk131 = my_rates->k13[logTlininterp_buf.indixe[i]-1] +
    (my_rates->k13[logTlininterp_buf.indixe[i]+1-1] -my_rates->k13[logTlininterp_buf.indixe[i]-1])*eqtdef;
  double heq1 = (-1. / (4.*eqk221)) * (eqk131-
    std::sqrt(8.*eqk131*eqk221*
              my_chemistry->HydrogenFractionByMass*local_rho+std::pow(eqk131,2.)));

  double dheq = (std::fabs(heq2-heq1)/(std::exp(eqt2) - std::exp(eqt1)))
    * (tgas[i]/p2d[i]) * edot[i];
  double heq = (-1. / (4.*k22[i])) * (k13[i]-
    std::sqrt(8.*k13[i]*k22[i]*
              my_chemistry->HydrogenFractionByMass*local_rho+std::pow(k13[i],2.)));

  return heq / dheq;
}

// -------------------------------------------------------------

/// Updates the iteration mask in the case where the user has specified that we
/// are using grackle as a part of a coupled radiative transfer calculation
///
/// @param[out] itmask the mask that will be overriden
/// @param[in] idx_range specifies the index-range
/// @param[in] kphHI view of the HI photo-ionization rate field
/// @param[in] my_chemistry specifies grackle settings (we probably don't need
///     to pass in everything)
static inline void coupled_rt_modify_itmask_(
  gr_mask_type* itmask,
  IndexRange idx_range,
  grackle::impl::View<gr_float***> kphHI,
  const chemistry_data* my_chemistry
)
{
  // adjust iteration mask if the caller indicates that they're using
  // Grackle in a coupled radiative-transfer/chemistry-energy calculation
  // (that has intermediate steps)
  if (my_chemistry->use_radiative_transfer == 1 &&
      my_chemistry->radiative_transfer_coupled_rate_solver == 1)  {
    // we only define behavior for radiative_transfer_intermediate_step
    // values of 0 or 1

    const int j = idx_range.j;
    const int k = idx_range.k;

    if (my_chemistry->radiative_transfer_intermediate_step == 1) {
      // the caller has invoked this chemistry-energy solver as an
      // intermediate step of a coupled radiative-transfer/chemistry-energy
      // calculation and they only want the solver consider cells where
      // the radiation is non-zero
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = (kphHI(i,j,k) > 0) ? MASK_TRUE : MASK_FALSE;
      }
    } else if (my_chemistry->radiative_transfer_intermediate_step == 0) {
      // the caller has invoked this chemistry-energy solver outside
      // of their coupled radiative-transfer/chemistry-energy calculation.
      // They want to apply the solver to cells where radiation is 0 (i.e.
      // locations where skipped by the coupled calculation)
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = (kphHI(i,j,k) > 0) ? MASK_FALSE : MASK_TRUE;
      }
    }
  }
}


// -------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int solve_rate_cool_g(
  int imetal, double dt, InternalGrUnits internalu,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, photo_rate_storage* my_uvb_rates
)
{

  // Density, energy and velocity fields fields

  grackle::impl::View<gr_float***> de(my_fields->e_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(my_fields->internal_energy, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Radiative transfer fields

  grackle::impl::View<gr_float***> kphHI(my_fields->RT_HI_ionization_rate, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Constants

#ifdef GRACKLE_FLOAT_4
  const gr_float tolerance = (gr_float)(1.0e-05);
#else
  const gr_float tolerance = (gr_float)(1.0e-10);
#endif

  const double mh_local_var = mh_grflt;
  const double pi_local_var = pi_fortran_val;

  // row temporaries

  std::vector<double> dtit(my_fields->grid_dimension[0]);
  std::vector<double> ttot(my_fields->grid_dimension[0]);
  std::vector<double> p2d(my_fields->grid_dimension[0]);
  std::vector<double> tgas(my_fields->grid_dimension[0]);
  std::vector<double> tdust(my_fields->grid_dimension[0]);
  std::vector<double> metallicity(my_fields->grid_dimension[0]);
  std::vector<double> dust2gas(my_fields->grid_dimension[0]);
  std::vector<double> rhoH(my_fields->grid_dimension[0]);
  std::vector<double> mmw(my_fields->grid_dimension[0]);
  std::vector<double> ddom(my_fields->grid_dimension[0]);

  // Rate equation row temporaries

  std::vector<double> dedot(my_fields->grid_dimension[0]);
  std::vector<double> HIdot(my_fields->grid_dimension[0]);
  std::vector<double> dedot_prev(my_fields->grid_dimension[0]);
  std::vector<double> HIdot_prev(my_fields->grid_dimension[0]);
  std::vector<double> k13dd(my_fields->grid_dimension[0] * 14);
  std::vector<double> h2dust(my_fields->grid_dimension[0]);

  // Cooling/heating row locals

  std::vector<double> edot(my_fields->grid_dimension[0]);

  // Iteration mask

  gr_mask_type anydust;
  std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);
  std::vector<gr_mask_type> itmask_tmp(my_fields->grid_dimension[0]);
  std::vector<gr_mask_type> itmask_nr(my_fields->grid_dimension[0]);
  std::vector<gr_mask_type> itmask_metal(my_fields->grid_dimension[0]);
  std::vector<int> imp_eng(my_fields->grid_dimension[0]);

  
  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  // Set error indicator (we will return this value)
  int ierr = GR_SUCCESS;

  // Set flag for dust-related options

  if ((my_chemistry->h2_on_dust > 0)  ||  (my_chemistry->dust_chemistry > 0))  {
    anydust = MASK_TRUE;
  } else {
    anydust = MASK_FALSE;
  }

  // ignore metal chemistry/cooling below this metallicity
  const double min_metallicity = 1.e-9 / my_chemistry->SolarMetalFractionByMass;
      
  // Set units
  const double dom      = internalu_calc_dom_(internalu);
  const double chunit   = internalu_get_chunit_(internalu);

  const double dx_cgs = my_fields->grid_dx * internalu.xbase1;
  const double c_ljeans = internalu_calc_coef_ljeans_(internalu,
                                                      my_chemistry->Gamma);

  const double dlogtem = (
    (std::log(my_chemistry->TemperatureEnd) -
     std::log(my_chemistry->TemperatureStart)) /
    (double)(my_chemistry->NumberOfTemperatureBins-1 )
  );

  // We better make consistent at first GC202002

  if (my_chemistry->primordial_chemistry > 0)  {

#define ABUNDANCE_CORRECTION
#ifdef ABUNDANCE_CORRECTION
    wrapped_make_consistent_g_(imetal, dom, my_chemistry, my_rates, my_fields);
#endif

  }

  // Convert densities from comoving to proper

  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)(std::pow(internalu.a_value,(-3)) );
    wrapped_scale_fields_g_(imetal, factor, my_chemistry, my_fields);
  }

#ifdef ABUNDANCE_CORRECTION
  wrapped_ceiling_species_g_(imetal, my_chemistry, my_fields);
#endif

  const grackle_index_helper idx_helper = build_index_helper_(my_fields);

  // parallelize the k and j loops with OpenMP
  // flat j and k loops for better parallelism
  //_// PORT: #ifdef _OPENMP
  //_// PORT: ! ierr is declared as shared and should be modified with atomic operation
  //_// PORT: !$omp parallel do schedule(runtime) private(
  //_// PORT: !$omp&   dtit, ttot, p2d, tgas,
  //_// PORT: !$omp&   tdust, metallicity, dust2gas, rhoH, mmw,
  //_// PORT: !$omp&   ddom,
  //_// PORT: !$omp&   dep, dedot,HIdot, dedot_prev,
  //_// PORT: !$omp&   HIdot_prev,
  //_// PORT: !$omp&   k13dd, h2dust,
  //_// PORT: !$omp&   edot,
  //_// PORT: !$omp&   itmask, itmask_metal )
  //_// PORT: #endif
  //_// TODO_USE: OMP_PRAGMA("omp parallel")
  {
    // TODO: move more relevant variable declarations to here to replace the
    //       OMP private-clause

    // holds computed grain temperatures:
    grackle::impl::GrainSpeciesCollection grain_temperatures =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    // holds computed grain growth/destruction rates:
    grackle::impl::GrainSpeciesCollection grain_growth_rates =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf =
      grackle::impl::new_LogTLinInterpScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf =
      grackle::impl::new_Cool1DMultiScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::CoolHeatScratchBuf coolingheating_buf =
      grackle::impl::new_CoolHeatScratchBuf(my_fields->grid_dimension[0]);

    // buffers in the following data structure are used to temporarily hold
    // the evolved density of various species as we evolve over a subcycle
    grackle::impl::SpeciesCollection species_tmpdens =
      grackle::impl::new_SpeciesCollection(my_fields->grid_dimension[0]);

    // buffers in the following data structure are used to temporarily hold
    // the interpolated Collisional/Recombination Rates that have interpolated
    // using the standard 1D log temperature table.
    grackle::impl::ColRecRxnRateCollection kcr_buf =
      grackle::impl::new_ColRecRxnRateCollection(my_fields->grid_dimension[0]);

    // buffers in the following data structure are used to temporarily hold
    // the computed radiative reaction rates
    grackle::impl::PhotoRxnRateCollection kshield_buf =
      grackle::impl::new_PhotoRxnRateCollection(my_fields->grid_dimension[0]);

    // buffers in the following data structure are used to temporarily hold
    // the interpolated chemistry-heating rates at each islice zone
    grackle::impl::ChemHeatingRates chemheatrates_buf =
      grackle::impl::new_ChemHeatingRates(my_fields->grid_dimension[0]);

    //_// TODO_USE: OMP_PRAGMA("omp for")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      // construct an index-range corresponding to "i-slice"
      const IndexRange idx_range = make_idx_range_(t, &idx_helper);
      const int k = idx_range.k; // use 0-based index
      const int j = idx_range.j; // use 0-based index

      // tolerance = 1.0e-06_DKIND * dt

      // Initialize iteration mask to true for all cells.

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = MASK_TRUE;
      }

      // adjust iteration mask (but only if using Grackle in a coupled
      // radiative-transfer calculation)
      coupled_rt_modify_itmask_(itmask.data(), idx_range, kphHI, my_chemistry);

      // Set time elapsed to zero for each cell in 1D section

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        ttot[i] = 0.;
      }

      // A useful slice variable since we do this a lot

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        ddom[i] = d(i,j,k) * dom;
      }

      // declare 2 variables (primarily used for subcycling, but also used in
      // error reporting)
      int iter;
      double ttmin;

      // ------------------ Loop over subcycles ----------------

      for (iter = 1; iter<=(my_chemistry->max_iterations); iter++) {

        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          if (itmask[i] != MASK_FALSE)  {
            dtit[i] = huge8;
          }
        }

        // Compute the cooling rate, tgas, tdust, and metallicity for this row

        wrapped_cool1d_multi_g_(
          imetal, idx_range, iter, edot.data(), tgas.data(),
          mmw.data(), p2d.data(), tdust.data(), metallicity.data(),
          dust2gas.data(), rhoH.data(), itmask.data(), itmask_metal.data(),
          my_chemistry, my_rates, my_fields, *my_uvb_rates, internalu,
          grain_temperatures, logTlininterp_buf, cool1dmulti_buf,
          coolingheating_buf
        );

        if (my_chemistry->primordial_chemistry == 0)  {
          // This is some basic book-keeping to ensure that itmask_tmp has
          // sensible values when ispecies is 0
          // -> There's room for optimization: when ispecies is 0, there is
          //    need to ever touch this variable. But, for now we focus on
          //    correct behavior before implementing this optimization
          std::memcpy(itmask_tmp.data(), itmask.data(), sizeof(gr_mask_type)*my_fields->grid_dimension[0]);

        } else {

          // Look-up rates as a function of temperature for 1D set of zones
          //  (maybe should add itmask to this call)

          wrapped_lookup_cool_rates1d_g_(
            idx_range, anydust, tgas.data(), mmw.data(),
            tdust.data(), dust2gas.data(), k13dd.data(), h2dust.data(),
            dom, dx_cgs, c_ljeans, itmask.data(), itmask_metal.data(),
            imetal, rhoH.data(), dt, my_chemistry, my_rates, my_fields,
            *my_uvb_rates, internalu, grain_growth_rates, grain_temperatures,
            logTlininterp_buf, kcr_buf, kshield_buf, chemheatrates_buf
          );

          // Compute dedot and HIdot, the rates of change of de and HI
          //   (should add itmask to this call)

          wrapped_rate_timestep_g_(
            dedot.data(), HIdot.data(), anydust, idx_range, h2dust.data(),
            rhoH.data(), itmask.data(), edot.data(), chunit, dom,
            my_chemistry, my_fields, *my_uvb_rates, kcr_buf, kshield_buf,
            chemheatrates_buf
          );

          // move itmask temporary array
          // then split cells with low densities
          //    => Gauss-Seidel scheme
          //              and with high densities
          //    => Newton-Raphson scheme

          std::memcpy(itmask_tmp.data(), itmask.data(), sizeof(gr_mask_type)*my_fields->grid_dimension[0]);
          std::memcpy(itmask_nr.data(), itmask.data(), sizeof(gr_mask_type)*my_fields->grid_dimension[0]);
          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            if ( itmask_tmp[i] != MASK_FALSE )  {

              if ( (imetal == 0) && (ddom[i] < 1.e8) ) {
                itmask_nr[i] = MASK_FALSE;
              } else if (
                (imetal == 1) &&
                (((metallicity[i] <= min_metallicity) && (ddom[i] < 1.e8)) ||
                 ((metallicity[i] > min_metallicity) && (ddom[i] < 1.e6)))
              ) {
                itmask_nr[i] = MASK_FALSE;
              } else {
                itmask[i] = MASK_FALSE;
              }

            }
          }

          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            if (itmask_nr[i] != MASK_FALSE)  {
              if ((my_chemistry->with_radiative_cooling == 1) &&
                  (my_chemistry->primordial_chemistry > 1) &&
                  (ddom[i] > 1.e7) && (tgas[i] > 1650.0)) {
                imp_eng[i] = 1;
              } else {
                imp_eng[i] = 0;
              }
            }
          }

          // Find timestep that keeps relative chemical changes below 10%

          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            if (itmask[i] != MASK_FALSE)  {
              // Bound from below to prevent numerical errors
               
              if (std::fabs(dedot[i]) < tiny8)
                         { dedot[i] = std::fmin(tiny_fortran_val,de(i,j,k)); }
              if (std::fabs(HIdot[i]) < tiny8)
                         { HIdot[i] = std::fmin(tiny_fortran_val,HI(i,j,k)); }

              // If the net rate is almost perfectly balanced then set
              //     it to zero (since it is zero to available precision)

              if (std::fmin(std::fabs(kcr_buf.data[ColRecRxnLUT::k1][i]* de(i,j,k)*HI(i,j,k)),
                      std::fabs(kcr_buf.data[ColRecRxnLUT::k2][i]*HII(i,j,k)*de(i,j,k)))/
                  std::fmax(std::fabs(dedot[i]),std::fabs(HIdot[i])) >
                   1.0e6)  {
                dedot[i] = tiny8;
                HIdot[i] = tiny8;
              }

              // If the iteration count is high then take the smaller of
              //   the calculated dedot and last time step's actual dedot.
              //   This is intended to get around the problem of a low
              //   electron or HI fraction which is in equilibrium with high
              //   individual terms (which all nearly cancel).

              if (iter > 50)  {
                dedot[i] = std::fmin(std::fabs(dedot[i]), std::fabs(dedot_prev[i]));
                HIdot[i] = std::fmin(std::fabs(HIdot[i]), std::fabs(HIdot_prev[i]));
              }

              // compute minimum rate timestep

              double olddtit = dtit[i];
              dtit[i] = grackle::impl::fmin(std::fabs(0.1*de(i,j,k)/dedot[i]),
                   std::fabs(0.1*HI(i,j,k)/HIdot[i]),
                   dt-ttot[i], 0.5*dt);

              if (ddom[i] > 1.e8  && 
                   edot[i] > 0.       && 
                  my_chemistry->primordial_chemistry > 1)  {
                // here, we ensure that that the equilibrium mass density of
                // Hydrogen changes by 10% or less
                double Heq_div_dHeqdt = calc_Heq_div_dHeqdt_(
                  my_chemistry, my_rates, dlogtem, logTlininterp_buf,
                  kcr_buf.data[ColRecRxnLUT::k13], kcr_buf.data[ColRecRxnLUT::k22], d(i,j,k), tgas.data(),
                  p2d.data(), edot.data(), i
                );

                dtit[i] = std::fmin(dtit[i], 0.1*Heq_div_dHeqdt);
              }
              if (iter>10LL)  {
                dtit[i] = std::fmin(olddtit*1.5, dtit[i]);
              }

            } else if ((itmask_nr[i]!=MASK_FALSE) && 
                     (imp_eng[i]==0))  {
              dtit[i] = grackle::impl::fmin(std::fabs(0.1*e(i,j,k)/edot[i]*d(i,j,k)),
                   dt-ttot[i], 0.5*dt);

            } else if ((itmask_nr[i]!=MASK_FALSE) && 
                     (imp_eng[i]==1))  {
              dtit[i] = dt - ttot[i];

            } else {
              dtit[i] = dt;
            }
          }

        }

        // Compute maximum timestep for cooling/heating

        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          if (itmask[i] != MASK_FALSE)  {
            // Set energy per unit volume of this cell based in the pressure
            // (the gamma used here is the right one even for H2 since p2d
            //  is calculated with this gamma).

            double energy = std::fmax(p2d[i]/(my_chemistry->Gamma-1.), tiny8);

            // If the temperature is at the bottom of the temperature look-up
            // table and edot < 0, then shut off the cooling.

            if (tgas[i] <= 1.01*my_chemistry->TemperatureStart  && 
                 edot[i] < 0.)
                 { edot[i] = tiny8; }
            if (std::fabs(edot[i]) < tiny8) { edot[i] = tiny8; }

            // Compute timestep for 10% change

            dtit[i] = grackle::impl::fmin((double)(std::fabs(0.1*
              energy/edot[i]) ),
              dt-ttot[i], dtit[i]);

            if (dtit[i] != dtit[i])  {
              OMP_PRAGMA_CRITICAL
              {
                eprintf("HUGE dtit ::  %g %g %g %g %g %g %g\n",
                        energy,
                        edot [ i ],
                        dtit [ i ],
                        dt,
                        ttot [ i ],
                        std::fabs ( 0.1 * energy / edot [ i ] ),
                        (double) ( std::fabs ( 0.1 * energy / edot [ i ] ) ));
              }
            }

          }
        }

        // Update total and gas energy

        if (my_chemistry->with_radiative_cooling == 1)  {
          for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
            if (itmask[i] != MASK_FALSE)  {
              e(i,j,k)  = e(i,j,k) +
                      (gr_float)(edot[i]/d(i,j,k)*dtit[i] );

            }
          }
        }

        if (my_chemistry->primordial_chemistry > 0)  {

          // Solve rate equations with one linearly implicit Gauss-Seidel
          // sweep of a backward Euler method (for all cells specified by
          // itmask)

          wrapped_step_rate_g_(
            dtit.data(), idx_range, anydust, h2dust.data(), rhoH.data(),
            dedot_prev.data(), HIdot_prev.data(), itmask.data(),
            itmask_metal.data(), imetal, my_chemistry, my_fields,
            *my_uvb_rates, grain_growth_rates, species_tmpdens, kcr_buf,
            kshield_buf
          );

          // Solve rate equations with one linearly implicit Gauss-Seidel
          // sweep of a backward Euler method (for all cells specified by
          // itmask_nr)

          wrapped_step_rate_newton_raphson_(
            imetal, idx_range, iter, dom, chunit, dx_cgs, c_ljeans,
            dtit.data(), p2d.data(), tgas.data(), tdust.data(),
            metallicity.data(), dust2gas.data(), rhoH.data(), mmw.data(),
            h2dust.data(), edot.data(), anydust, itmask_nr.data(),
            itmask_metal.data(), imp_eng.data(), my_chemistry, my_rates,
            my_fields, *my_uvb_rates, internalu, grain_temperatures,
            logTlininterp_buf, cool1dmulti_buf, coolingheating_buf,
            chemheatrates_buf
          );

        }

        // return itmask
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          itmask[i] = itmask_tmp[i];
        }

        // Add the timestep to the elapsed time for each cell and find
        //  minimum elapsed time step in this row

        ttmin = huge8;
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          ttot[i] = std::fmin(ttot[i] + dtit[i], dt);
          if (std::fabs(dt-ttot[i]) <
               tolerance*dt) { itmask[i] = MASK_FALSE; }
          if (ttot[i]<ttmin) { ttmin = ttot[i]; }
        }

        // If all cells are done (on this slice), break out of subcycle loop
        if (std::fabs(dt-ttmin) < tolerance*dt) { break; }

      }  // subcycle iteration loop (for current idx_range)

      // review number of iterations that were spent in the subcycle loop

      if (iter > my_chemistry->max_iterations)  {
        OMP_PRAGMA_CRITICAL
        {
          printf("inside if statement solve rate cool: %d %d\n",
                 my_fields->grid_start[0],
                 my_fields->grid_end[0]);
          eprintf("MULTI_COOL iter >  %d  at j_0based,k_0based = %d %d\n",
                  my_chemistry->max_iterations, idx_range.j, idx_range.k);
          printf("FATAL error (2) in MULTI_COOL\n");
          printf("( dt = %.17e ttmin = %.17e )", dt, ttmin);
          grackle::impl::print_contiguous_row_(
            dtit.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            ttot.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            edot.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );
          grackle::impl::print_contiguous_row_(
            itmask.data(), my_fields->grid_start[0], my_fields->grid_end[0]+1
          );

          if (my_chemistry->exit_after_iterations_exceeded == 1)  {
            ierr = GR_FAIL;
          }
        }  // OMP_PRAGMA_CRITICAL
      }

      if (iter > my_chemistry->max_iterations/2) { // WARNING_MESSAGE
        OMP_PRAGMA_CRITICAL
        {
          eprintf("MULTI_COOL iter,j_0based,k_0based = %d %d %d\n",
                  iter, idx_range.j, idx_range.k);
        }
      }

    }  // outer-loop (index t) - each of these correspond to j,k pairs

    // cleanup manually allocated temporaries
    grackle::impl::drop_GrainSpeciesCollection(&grain_temperatures);
    grackle::impl::drop_GrainSpeciesCollection(&grain_growth_rates);
    grackle::impl::drop_LogTLinInterpScratchBuf(&logTlininterp_buf);
    grackle::impl::drop_Cool1DMultiScratchBuf(&cool1dmulti_buf);
    grackle::impl::drop_CoolHeatScratchBuf(&coolingheating_buf);
    grackle::impl::drop_SpeciesCollection(&species_tmpdens);
    grackle::impl::drop_ColRecRxnRateCollection(&kcr_buf);
    grackle::impl::drop_PhotoRxnRateCollection(&kshield_buf);
    grackle::impl::drop_ChemHeatingRates(&chemheatrates_buf);


  }  // OMP_PRAGMA("omp parallel")

  // If an error has been produced, return now.

  if (ierr != GR_SUCCESS)  {
    return ierr;
  }

  // Convert densities back to comoving from proper

  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)(std::pow(internalu.a_value,3) );
    wrapped_scale_fields_g_(imetal, factor, my_chemistry, my_fields);
  }

  if (my_chemistry->primordial_chemistry > 0)  {

    // Correct the species to ensure consistency (i.e. type conservation)

#ifdef ABUNDANCE_CORRECTION
    wrapped_make_consistent_g_(imetal, dom, my_chemistry, my_rates, my_fields);
    wrapped_ceiling_species_g_(imetal, my_chemistry, my_fields);
#endif

  }

  return ierr;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
