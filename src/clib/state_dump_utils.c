/***********************************************************************
/
/ Utilities related to dumping state (for debugging purposes)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#include <stdio.h>

#include "grackle.h"

// define some functionality to help write json objects
// ====================================================

struct json_obj_writer {
  FILE* fp;
  const char* member_spacing; // holds the whitespace between members
  int num_separations;
};

static struct json_obj_writer json_create_writer_(FILE *fp){
  fputc('{', fp); // intentionally omit '\n'
  struct json_obj_writer out = {fp, "\n  ", 0};
  return out;
}

void json_write_separation_(struct json_obj_writer* writer){
  if (writer->num_separations > 0) fputc(',', writer->fp);
  fprintf(writer->fp, "%s", writer->member_spacing);
  writer->num_separations++;
}

void json_finish_(struct json_obj_writer* writer){
  if (writer->num_separations == 0) {
    fprintf(writer->fp, "}\n");
  } else {
    fprintf(writer->fp, "\n}\n");
  }
}

// Define helpers for the show_parameters function
static void json_field_INT(struct json_obj_writer* writer, const char* field,
                           int val)
{
  json_write_separation_(writer);
  fprintf(writer->fp, "\"%s\" : %d", field, val);
}
static void json_field_DOUBLE(struct json_obj_writer* writer,
                              const char* field, double val)
{
  json_write_separation_(writer);
  fprintf(writer->fp, "\"%s\" : %.17g", field, val);
}
static void json_field_STRING(struct json_obj_writer* writer,
                              const char* field, const char* val)
{
  json_write_separation_(writer);
  if (val == NULL){
    fprintf(writer->fp, "\"%s\" : null", field);
  } else {
    fprintf(writer->fp, "\"%s\" : \"%s\"", field, val);
  }
}

// functions to write json representations of various structs to disk
// ==================================================================

void show_parameters_(FILE *fp, const chemistry_data *my_chemistry){
  struct json_obj_writer writer = json_create_writer_(fp);

  // write all of the fields
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) \
    json_field_ ## TYPE (&writer, #FIELD, my_chemistry->FIELD);
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY

  // end the json object
  json_finish_(&writer);
}

void show_version_(FILE *fp, const grackle_version* gversion)
{
  struct json_obj_writer writer = json_create_writer_(fp);
  json_field_STRING(&writer, "version", gversion->version);
  json_field_STRING(&writer, "branch", gversion->branch);
  json_field_STRING(&writer, "revision", gversion->revision);
  json_finish_(&writer);
}
