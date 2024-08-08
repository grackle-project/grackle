/***********************************************************************
/
/ Declare functions used to run the show-* subcommands.
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef GRCLI_CMDSHOW_H
#define GRCLI_CMDSHOW_H

#include "CliParser.h"

namespace cmd::show {

enum class Kind {
  initial_units,
  parameters,
  precision,
  version,
};

/// execute the show-* command
[[noreturn]] void run(CliParser& parser, Kind kind);

} // namespace cmd::show

#endif /* GRCLI_CMDSHOW_H */
