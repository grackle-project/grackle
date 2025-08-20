/// @file status_reporting.c
/// @brief Implements status reporting functionality

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "status_reporting.h"
#include "grackle.h" // GR_FAIL

// this is the internal routine that everything else dispatches to
static void vprint_err_(
  int internal_error, const struct grimpl_source_location_ locinfo,
  const char* msg, va_list vlist
) {
  const char* santized_func_name = (locinfo.fn_name == NULL)
    ? "{unspecified}" : locinfo.fn_name;

  const char* fallback_msg_ = "{NULL encountered instead of error message}";
  char* dynamic_msg_buf = NULL;
  const char* msg_buf;
  if (msg == NULL) {
    msg_buf = fallback_msg_;
  } else {
    // make a copy of the variadic function arguments
    va_list vlist_copy;
    va_copy(vlist_copy, vlist);

    // get the total size of the formatted message
    size_t msg_len = vsnprintf(NULL, 0, msg, vlist_copy) + 1;
    va_end(vlist_copy);

    // allocate the buffer to hold the message
    dynamic_msg_buf = malloc(sizeof(char) * msg_len);

    // actually format the message
    vsnprintf(dynamic_msg_buf, msg_len, msg, vlist);
    va_end(vlist);

    msg_buf = dynamic_msg_buf;
  }

  const char* descr = (internal_error == 1) ? "FATAL" : "ERROR";

  fprintf(
    stderr, "Grackle-%s %s:%d in %s] %s\n", descr, locinfo.file,
    locinfo.lineno, santized_func_name, msg_buf
  );

  if (dynamic_msg_buf != NULL) { free(dynamic_msg_buf); }
}

void grimpl_abort_with_internal_err_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
) {
  va_list args;
  va_start(args, msg);
  vprint_err_(1, locinfo, msg, args);
  // while stderr should flush by default, people may overwrite this behavior.
  // We want to force flushing here since we are aborting the program
  fflush(stderr);
  abort();
}

int grimpl_print_and_return_err_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
) {
  va_list args;
  va_start(args, msg);
  vprint_err_(0, locinfo, msg, args);
  return GR_FAIL;
}

void grimpl_print_err_msg_(
  const struct grimpl_source_location_ locinfo, const char* msg, ...
) {
  va_list args;
  va_start(args, msg);
  vprint_err_(0, locinfo, msg, args);
}

