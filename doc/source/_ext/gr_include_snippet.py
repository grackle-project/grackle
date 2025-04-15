"""
This is a simple extension, gr_include_snippet, that introduces a Directive
that acts like Sphinx's literalinclude directive

In more detail, the idea is to:
- have an external example that has specially-escaped anchors
- the gr-include-snippet is used to skip over these anchors
- we use the anchors to extract snippets of the code in literalinclude snippets
  to annotate parts of the example.
"""

from __future__ import annotations
import re

from docutils import nodes
from docutils.parsers.rst import directives

from sphinx.application import Sphinx
from sphinx.util.docutils import SphinxDirective

try:
    from sphinx.util.typing import ExtensionMetadata
except ImportError:
    from typing import Any

    ExtensionMetadata = Any

_IGNORE_REGEX = re.compile(r"^\s*" + re.escape(r"//@%"))


def _filtered_lines(fname, *, encoding=None):
    # produces an iterator over valid lines in fname (trailing whitespace is removed)
    yielded_any = False
    with open(fname, encoding=encoding) as f:
        for line in f:
            if _IGNORE_REGEX.match(line) is not None:
                continue
            elif line.isspace() and not yielded_any:  # skip any leading whitespace
                continue
            yielded_any = True
            yield line.rstrip()


class GrIncludeSnippetDirective(SphinxDirective):
    required_arguments = 1
    optional_arguments = 0
    option_spec = {
        "encoding": directives.encoding,
        "language": directives.unchanged_required,
        # below, we use directives.unchanged to make forwarding easier
        "caption": directives.unchanged,
        "class": directives.unchanged,
        "name": directives.unchanged,
    }

    def run(self) -> list[nodes.Node]:
        # get path to the file that will be read
        rel_path_to_docroot, abs_path = self.env.relfn2path(self.arguments[0])
        self.env.note_dependency(rel_path_to_docroot)

        lines = []
        _indent = "   "
        lines.append(f".. code-block:: {self.options['language']}")
        for optname in self.option_spec:
            val = self.options.get(optname, None)
            if optname in ["encoding", "language"] or val is None:
                continue
            elif val == "":
                lines.append(f"{_indent}:{optname}:")
            else:
                lines.append(f"{_indent}:{optname}: {val}")
        lines.append(f"{_indent}")

        it = _filtered_lines(abs_path, encoding=self.options.get("encoding", None))
        for line in it:
            lines.append(f"{_indent}{line}")

        # the following line is inspired by yt:
        self.state_machine.insert_input(
            input_lines=lines, source=self.state_machine.document.attributes["source"]
        )
        return []


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_directive("gr-include-snippet", GrIncludeSnippetDirective)

    return {"version": "0.1", "parallel_read_safe": True, "parallel_write_save": False}
