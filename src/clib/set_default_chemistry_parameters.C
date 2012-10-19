#include <string.h>
#include <stdio.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"

chemistry_data set_default_chemistry_parameters()
{
  
  chemistry_data my_chemistry;

  my_chemistry.Gamma                          = 5./3.;
  my_chemistry.use_chemistry                  = FALSE;  // off
  my_chemistry.primordial_chemistry           = FALSE;  // off
  my_chemistry.metal_cooling                  = FALSE;
  my_chemistry.h2_on_dust                     = FALSE;

  my_chemistry.cmb_temperature_floor          = TRUE;   // use CMB floor.
  my_chemistry.cloudy_table_file              = "";
  my_chemistry.include_metal_heating          = FALSE;

  my_chemistry.three_body_rate                = 0;   // ABN02
  my_chemistry.cie_cooling                    = 1;
  my_chemistry.h2_optical_depth_approximation = 1;

  my_chemistry.photoelectric_heating          = 0;
  my_chemistry.photoelectric_heating_rate     = 8.5e-26;  // ergs cm-3 s-1

  my_chemistry.RadiationFieldType             = 0;
  my_chemistry.RadiationFieldRedshift         = 0.0;
  my_chemistry.TabulatedLWBackground          = 0;
  my_chemistry.RadiationFieldLevelRecompute   = 0;
  my_chemistry.AdjustUVBackground             = 1;
  my_chemistry.AdjustUVBackgroundHighRedshift = 0;
  my_chemistry.SetUVBAmplitude                = 1.0;
  my_chemistry.SetHeIIHeatingScale            = 1.8;
  my_chemistry.RadiationXRaySecondaryIon      = 0;
  my_chemistry.RadiationXRayComptonHeating    = 0;

  my_chemistry.f3                       = 1.0e-21;        // radiation normalization
  my_chemistry.f0to3                    = 0.1;
  my_chemistry.RadiationRedshiftOn      = 7.0;
  my_chemistry.RadiationRedshiftOff     = 0.0;
  my_chemistry.RadiationRedshiftFullOn  = 6.0;
  my_chemistry.RadiationRedshiftDropOff = 0.0;
  my_chemistry.HydrogenFractionByMass   = 0.76;
  /* The DToHRatio is by mass in the code, so multiply by 2. */
  my_chemistry.DeuteriumToHydrogenRatio = 2.0*3.4e-5; // Burles & Tytler 1998
  my_chemistry.SolarMetalFractionByMass = 0.02041;
  my_chemistry.NumberOfTemperatureBins = 600;
  my_chemistry.ih2co                   = 1;
  my_chemistry.ipiht                   = 1;
  my_chemistry.TemperatureStart        = 1.0;
  my_chemistry.TemperatureEnd          = 1.0e8;
  my_chemistry.comp_xray               = 0;
  my_chemistry.temp_xray               = 0;
  my_chemistry.CaseBRecombination      = 0;   // default to case A rates
  my_chemistry.NumberOfDustTemperatureBins = 250;
  my_chemistry.DustTemperatureStart    = 1.0;
  my_chemistry.DustTemperatureEnd      = 1500.0;

  my_chemistry.CloudyCoolingGridRank          = 0;
  my_chemistry.CloudyElectronFractionFactor = 9.153959e-3; // Cloudy 07.02 abundances

  return my_chemistry;
}
