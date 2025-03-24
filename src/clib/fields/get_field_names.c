// See LICENSE file for license and copyright information

/// @file get_field_names.c
/// @brief Defines the function to get the list of field names.

#include "FrozenKeyIdxBiMap.h"
#include "get_field_names.h"
#include "grackle.h"
#include "status_reporting.h"

enum FieldReq {
  FIELDREQ_NEVER_SATISFIED,
  FIELDREQ_ALWAYS_SATISFIED,
  FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1,
  FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_2,
  FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_3,
  FIELDREQ_METAL_COOLING_EQ_1,
  FIELDREQ_DUST_CHEMISTRY_EQ_1,
  FIELDREQ_USE_VOLUMETRIC_HEATING_RATE,
  FIELDREQ_USE_SPECIFIC_HEATING_RATE,
  FIELDREQ_USE_TEMPERATURE_FLOOR,
  FIELDREQ_USE_RADIATIVE_TRANSFER,
  FIELDREQ_H2_SELF_SHIELDING_EQ_2,
  FIELDREQ_H2_CUSTOM_SHIELDING_EQ_1,
  FIELDREQ_USE_ISRF_FIELD,
};

// returns 1 for true, 0 for false, -1 for unrecognized req
static int meets_req_(const chemistry_data* my_chemistry, enum FieldReq req)
{
  switch (req) {
  case FIELDREQ_NEVER_SATISFIED:
    return 0;
  case FIELDREQ_ALWAYS_SATISFIED:
    return 1;
  case FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1:
    return (my_chemistry->primordial_chemistry >= 1);
  case FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_2:
    return (my_chemistry->primordial_chemistry >= 2);
  case FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_3:
    return (my_chemistry->primordial_chemistry >= 3);
  case FIELDREQ_METAL_COOLING_EQ_1:
    return (my_chemistry->metal_cooling == 1);
  case FIELDREQ_DUST_CHEMISTRY_EQ_1:
    return (my_chemistry->dust_chemistry == 1);
  case FIELDREQ_USE_VOLUMETRIC_HEATING_RATE:
    return (my_chemistry->use_volumetric_heating_rate == 1);
  case FIELDREQ_USE_SPECIFIC_HEATING_RATE:
    return (my_chemistry->use_specific_heating_rate == 1);
  case FIELDREQ_USE_TEMPERATURE_FLOOR:
    return (my_chemistry->use_temperature_floor == 1);
  case FIELDREQ_USE_RADIATIVE_TRANSFER:
    return (my_chemistry->use_radiative_transfer == 1);
  case FIELDREQ_H2_SELF_SHIELDING_EQ_2:
    return (my_chemistry->H2_self_shielding == 2);
  case FIELDREQ_H2_CUSTOM_SHIELDING_EQ_1:
    return (my_chemistry->H2_custom_shielding == 1);
  case FIELDREQ_USE_ISRF_FIELD:
    return (my_chemistry->use_isrf_field == 1);
  default:
    return -1;
  }
}

struct field_name_req_pair_ {
  const char* name;
  enum FieldReq req;
};

static const struct field_name_req_pair_ field_req_l_[] = {
  // in the future, we could probably use XMacros to build this
  {"density", FIELDREQ_ALWAYS_SATISFIED},
  {"internal_energy", FIELDREQ_ALWAYS_SATISFIED},
  {"x_velocity", FIELDREQ_NEVER_SATISFIED},
  {"y_velocity", FIELDREQ_NEVER_SATISFIED},
  {"z_velocity", FIELDREQ_NEVER_SATISFIED},

  {"HI_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1},
  {"HII_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1},
  {"HeI_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1},
  {"HeII_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1},
  {"HeIII_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1},
  {"e_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_1},

  {"H2I_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_2},
  {"H2II_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_2},
  {"HM_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_2},

  {"DI_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_3},
  {"DII_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_3},
  {"HDI_density", FIELDREQ_PRIMORDIAL_CHEMISTRY_GE_3},

  {"metal_density", FIELDREQ_METAL_COOLING_EQ_1},
  {"dust_density", FIELDREQ_DUST_CHEMISTRY_EQ_1},

  {"volumetric_heating_rate", FIELDREQ_USE_VOLUMETRIC_HEATING_RATE},
  {"specific_heating_rate", FIELDREQ_USE_SPECIFIC_HEATING_RATE},
  {"temperature_floor", FIELDREQ_USE_TEMPERATURE_FLOOR},

  {"RT_heating_rate", FIELDREQ_USE_RADIATIVE_TRANSFER},
  {"RT_HI_ionization_rate", FIELDREQ_USE_RADIATIVE_TRANSFER},
  {"RT_HeI_ionization_rate", FIELDREQ_USE_RADIATIVE_TRANSFER},
  {"RT_HeII_ionization_rate", FIELDREQ_USE_RADIATIVE_TRANSFER},
  {"RT_H2_dissociation_rate", FIELDREQ_USE_RADIATIVE_TRANSFER},

  {"H2_self_shielding_length", FIELDREQ_H2_SELF_SHIELDING_EQ_2},
  {"H2_custom_shielding_factor", FIELDREQ_H2_CUSTOM_SHIELDING_EQ_1},
  {"isrf_habing", FIELDREQ_USE_ISRF_FIELD}
};

static const int field_req_l_len_ = (
  sizeof(field_req_l_) / sizeof(struct field_name_req_pair_)
);


int get_field_names_from_local_(FrozenKeyIdxBiMap** out,
                                const chemistry_data* my_chemistry)
{
  // graceful arg checking should happen further up the call-stack (this is
  // just a sanity check!)
  GR_INTERNAL_REQUIRE((out != NULL) && (my_chemistry != NULL),
                      "Argument sanity checks failed");

  const char* field_name_l[field_req_l_len_];
  int name_l_count = 0;

  for (int i = 0; i < field_req_l_len_; i++) {
    int req_is_satisfied = meets_req_(my_chemistry, field_req_l_[i].req);
    GR_INTERNAL_REQUIRE(
      req_is_satisfied != -1,
      "the meets_req_ function doesn't seem to understand the requirement "
      "associated with the `%s` field",
      field_req_l_[i].name
    );

    if (req_is_satisfied) {
      field_name_l[name_l_count] = field_req_l_[i].name;
      name_l_count++;
    }
  }

  GR_INTERNAL_REQUIRE(
    name_l_count >= 1,
    "It shouldn't be possible for name_l_count to be non-positive"
  );

  // we use BIMAP_REFS_KEYDATA since all of the names are string literals
  return new_FrozenKeyIdxBiMap(out, field_name_l, name_l_count,
                               BIMAP_REFS_KEYDATA);
}
