#include <string.h>
#include <stdio.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"

int set_default_chemistry_data(chemistry_data &my_chemistry)
{
  
  my_chemistry.RadiativeCooling            = FALSE;             // off
  my_chemistry.MultiSpecies                = FALSE;             // off
  my_chemistry.ThreeBodyRate               = 0;                 // ABN02
  my_chemistry.CIECooling                  = 1;
  my_chemistry.H2OpticalDepthApproximation = 1;
  my_chemistry.H2FormationOnDust           = FALSE;


  my_chemistry.RadiationFieldType          = 0;
  my_chemistry.RadiationFieldRedshift      = 0.0;
  my_chemistry.TabulatedLWBackground       = 0;
  my_chemistry.RadiationFieldLevelRecompute = 0;
  my_chemistry.RadiationShield = 0;
  my_chemistry.AdjustUVBackground          = 1;
  my_chemistry.AdjustUVBackgroundHighRedshift = 0;
  my_chemistry.SetUVBAmplitude             = 1.0;
  my_chemistry.SetHeIIHeatingScale         = 1.8;
  my_chemistry.PhotoelectricHeating	      = 0;
  my_chemistry.PhotoelectricHeatingRate    = 8.5e-26;           // ergs cm-3 s-1
  my_chemistry.RadiationXRaySecondaryIon   = 0;
  my_chemistry.RadiationXRayComptonHeating = 0;

  my_chemistry.alpha0             = 1.5;               // radiation spectral slope
  my_chemistry.f3                 = 1.0e-21;           // radiation normalization
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
  my_chemistry.CloudyCoolingGridFile          = "";
  my_chemistry.IncludeCloudyHeating           = 0;
  my_chemistry.CMBTemperatureFloor            = 1;         // use CMB floor.
  my_chemistry.CloudyElectronFractionFactor = 9.153959e-3; // calculated using Cloudy 07.02 abundances

  my_chemistry.MetalCooling = FALSE;
  my_chemistry.MetalCoolingTable = (char*) "metal_cool.dat";

  return SUCCESS;
}
