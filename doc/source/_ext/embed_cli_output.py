"""
A simple sphinx extension that embeds output of a command line program
"""

from typing import Any, ClassVar, List, TYPE_CHECKING

import os
import subprocess
import sys

from docutils import nodes
from docutils.parsers.rst import directives

from sphinx.util.docutils import SphinxDirective

if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.util.typing import ExtensionMetadata, OptionSpec
else:
    Sphinx = Any
    ExtensionMetadata = Any
    OptionSpec = Any


class EmbedCLIOutput(SphinxDirective):
    """
    Implements ``.. embed-cli-output:: <path/to/prog>``

    You can use ``:args:`` to specify the arguments that are passed to the program
    ::
      .. embed-cli-output:: <path/to/prog>
         :args: --argA A --arg2

    The ``:literal:`` flag provides embeds the output as a literal block
    ::
      .. embed-cli-output:: <path/to/prog>
         :literal:
    """

    has_content = False
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False

    option_spec: ClassVar[OptionSpec] = {
        "args": directives.unchanged_required,
        "literal": directives.flag,
    }

    def run(self) -> List[nodes.Node]:
        # get absolute path to program
        _, path = self.env.relfn2path(self.arguments[0].strip())
        if os.path.isfile(path):
            self.state.document.settings.record_dependencies.add(path)
        else:
            self.severe(f"can't locate ``{path}``")

        # prepare the command
        cmd = [sys.executable, str(path)] if path.endswith(".py") else [str(path)]
        if self.options.get("args"):  # if diff is set, set udiff
            args = self.options["args"].strip().split()
            if len(args) == 0:
                self.severe("sanity check failed: `:args:` specified without args")
            else:
                cmd.extend(args)

        # run the command and prepare output
        rslt = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if rslt.returncode == 0:
            content = rslt.stdout.decode()
            force_literal_block = False
        else:
            msg = f"there was a problem executing ``{' '.join(cmd)}``"
            self.document.reporter.warning(msg, line=self.lineno)
            content = f"<{msg}>"
            force_literal_block = True
        content = content.strip()

        # return the output
        if force_literal_block or "literal" in self.options:
            return [nodes.literal_block(content)]
        else:
            return self.parse_text_to_nodes(
                content,
                allow_section_headings=False,  # <- this is the default
            )


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_directive("embed-cli-output", EmbedCLIOutput)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
