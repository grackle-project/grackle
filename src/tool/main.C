#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string_view>

#include "CliParser.h"
#include "cmd_bench.h"
#include "cmd_show.h"
#include "utils.h"


static void print_usage_(FILE* f, const char* bin_name) {
  std::fprintf(f,
      R"HELP(
USAGE:

    %s bench -S... PARAMETER_OPTIONS
    %s show-version
    %s show-precision
    %s show-parameters PARAMETER_OPTIONS
)HELP",
      bin_name, bin_name, bin_name, bin_name);
}

static const char* help_message_ = R"HELP(
DESCRIPTION:
    A simple command-line tool intended for debugging and basic benchmarking of
    grackle.


    The different modes are currently experimental:
      -> bench sets up an example problem and tries to run a benchmark
      -> show-parameters shows the parameters... (this may change a little in the future)

OPTIONS:
    Some Flags are grouped into categories

    -S
        Prefix for all "scenario arguments".

        At the moment, this might include
          * -Sgrid.ax<I>=<quantity>;<seq>
        Where the quantities in brackets get replaced by various things. For
        example:
          * <I> is replaced with 0, 1, or 2
          * <quantity> is replaced with temperature, density, metallicity
          * <seq> is replaced by arguments that specify a sequence of 1 or more
            values. In practice, this is either a lone floating-point value OR
            geomspace;<start>;<end>;<num>

        
    -h; --help
        Prints this message

  PARAMETER OPTIONS

    --par-start
        Every "word" after this option (until the --par-stop option) is parsed
        as a key-value pair that is used as a grackle parameter.

    --par-stop
        The sentinel value that must follow --par-start
)HELP";

/// prints the help message
///
/// @note
/// when somebody passes an -h flag, a help message is the expected result, so
/// we print to stdout (rather than stderr)
void print_help(const char* bin_name) {
  print_usage_(stdout, bin_name);
  std::printf("%s", help_message_);
}

int main(int argc, char* argv[]) {
  CliParser parser(argc, argv);

  // parse the first argument (i.e. the subcommand)
  const char* arg = parser.next();
  if (arg == nullptr) {
    std::fprintf(stderr, "%s is missing a subcommand", parser.bin_name());
    std::exit(1);
  }
  // handle the case where arg is a -h or --help
  try_parse_help(arg, parser.bin_name());

  std::string_view subcommand(arg);

  if (subcommand == "bench") {
    cmd::bench::run(parser);
  } else if (subcommand == "show-parameters") {
    cmd::show::run(parser, cmd::show::Kind::parameters);
  } else if (subcommand == "show-precision") {
    cmd::show::run(parser, cmd::show::Kind::precision);
  } else if (subcommand == "show-version") {
    cmd::show::run(parser, cmd::show::Kind::version);
  } else {
    std::fprintf(stderr, "unknown subcommand: %s\n", arg);
    std::exit(1);
  }

  GRCLI_ERROR("Something is wrong, this shouldn't be executed!");
}



