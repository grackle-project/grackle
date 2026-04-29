"""
Define logic pertaining to ConfReader
"""

import copy
from functools import reduce
import operator
import os
import string
import sys
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Collection,
    Dict,
    List,
    Mapping,
    NamedTuple,
    Sequence,
    Optional,
    Tuple,
    Type,
    Union,
)

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomli as tomllib
    except ImportError:
        import vendored_tomli as tomllib

if TYPE_CHECKING:
    import datetime

    _ScalarTOMLVal = Union[
        bool, str, float, int, datetime.datetime, datetime.time, datetime.date
    ]
    TOMLVal = Union[_ScalarTOMLVal, List["TOMLVal"], Dict[str, "TOMLVal"]]

else:
    TOMLVal = Any

# represents a string, specifying a key in the root table of a TOML document (i.e. the
# top level table), or a sequence of strings where each element is a key of a more
# nested table (directly analogous to a dotted TOML key).
TOMLPath = Union[str, Sequence[Union[str, int]]]


# we are going to be fairly pedantic about toml parsing


class ConfReadError(ValueError):
    """Subclass of ValueError used to provide a detailed error."""

    def __init__(self, msg: str, reader: "ConfReader", path: Optional[TOMLPath] = None):
        if path is None:
            err_loc = f"{reader.file_path}"
        else:
            path = (path,) if isinstance(path, str) else path
            parts = []
            for i, part in enumerate(path):
                if isinstance(part, int):
                    parts.append(f"[{i}]")
                else:
                    if i > 0:
                        parts.append(".")
                    parts.append(part)
            path = "".join(parts)
            err_loc = f"{path} in {reader.file_path}"
        ValueError.__init__(self, f"{msg}: {err_loc}")


_CheckTypeArg = Union[Type, Tuple[Type, ...], Callable[[Any], bool]]


def _check_type(v: Any, classinfo: _CheckTypeArg) -> bool:
    if isinstance(classinfo, (type, tuple)):
        return isinstance(v, classinfo)
    return classinfo(v)


class TOMLValReq(NamedTuple):
    """
    Specifies requirements for a toml value

    Checks are applied after the values are coerced to python values.
    """

    # this the table key
    name: str

    # this is either an argument for isinstance or a callable specifying whether
    # the value has the appropriate type
    classinfo: _CheckTypeArg

    # describes the required toml type
    type_descr: str

    description: str

    # when specified, the value is required to be one of these values
    choices: Optional[Collection] = None

    # whether a key is required
    required: bool = True

    def enforce(self, val: TOMLVal, reader: "ConfReader", path: TOMLPath):
        """Raise an error if ``val`` satisfy requirements."""

        if not _check_type(val, self.classinfo):
            msg = "invalid type, must be {self.type_descr}"
            raise ConfReadError(msg, reader=reader, path=path)
        elif (self.choices is not None) and (val not in self.choices):
            pretty_choices = ", ".join(repr(c) for c in self.choices)
            msg = f"invalid choice {val!r}, choose from {pretty_choices}"
            raise ConfReadError(msg, reader=reader, path=path)


# the following object primarily exists for exposition purposes:
# - ConfReader.process_table checks if its a specified req
# - ConfReader.read_parametrize then manually parses this parameter
PARAMETRIZED_TOML_REQ = TOMLValReq(
    name="parametrize",
    classinfo=list,
    type_descr="array of tables",
    required=False,
    description="""\
This parameter is used to actually parameterize a test case.
In more detail, this parameter is associated with an array of tables where each contained table associates the variable-name keys with values.
The variable-name keys must be consistent between tables.

For concreteness, lets consider the following example:

.. code-block:: toml

   parametrize = [
     {image="shared_build", overide_lib_path=false, extra_conf_arg=""},
     {image="static_build", overide_lib_path=false, extra_conf_arg="-GNinja"}
   ]

.. note::

   Reminder: the openning ``{`` and closing ``}`` of an inline-table must be on the
   same line in TOML 1.0. (This is relaxed in TOML 1.1)

The defined parametrized variables can be used to define other parameters (these other parameters must be sibling keys of parametrize) in 2 ways.

1. Directly reuse the value. This is illustrated in the following snippet

   .. code-block:: toml
   
      image.use-param = "image"
      override_LD_LIBRARY_PATH.use-param = "overide_lib_path"

2. String interpolation.
   We check all parsed strings for the python's
   `curly-brace format string syntax <https://docs.python.org/3/library/string.html#format-string-syntax>`__
   that is understood by ``str.format()``.
   The following snippet illustrates what this might look like

   .. code-block:: toml
   
      image = "{image}"
      cmds = ["cmake -Bbuild {extra_conf_arg}", "cmake --build build"]

   Reminder: ``{{`` escapes ``{`` and ``}`` escapes ``}``.""",
)


def _is_single_type_toml_arr(v: Any, kind: _CheckTypeArg, allow_empty: bool) -> bool:
    """Check if ``v`` corresponds to a toml array of a single type"""
    if not isinstance(v, list):
        return False
    elif len(v) == 0:
        return allow_empty
    return all(_check_type(elem, kind) for elem in v)


class ParametrizedProxy(NamedTuple):
    """A proxy for a parametrized value.

    This should be constructed with check_parametrized
    """

    # the nominal parameter being read
    primary_path: Tuple[str, ...]
    # keys are paths relative to primary_map where values get modified. The associated
    # value is a string if we perform a direct substitution or None in places where
    # we call .format(...) on the user provided strings
    sub_path_map: Dict[Tuple[Union[str, int], ...], Optional[str]]
    # the configuration requirement
    req: TOMLValReq

    def substitute(
        self, reader: "ConfReader", param_vars: Mapping[str, Any]
    ) -> TOMLVal:
        """Returns the corresponding value after substituting parametrized vals"""
        # make sure this method functions in a manner consistent with the description
        # provided PARAMETRIZED_TOML_REQ.description
        obj = copy.deepcopy(reader[self.primary_path])
        is_scalar = () in self.sub_path_map
        if is_scalar:
            obj = [obj]

        for rel_path, direct_sub_varname in self.sub_path_map.items():
            # print(rel_path, direct_sub_varname)
            path = self.primary_path if is_scalar else self.primary_path + rel_path
            if is_scalar:
                innermost_container = obj
                _loc = 0
            else:
                innermost_container = reduce(operator.getitem, rel_path[:-1], obj)
                _loc = rel_path[-1]
            raw_val = innermost_container[_loc]

            try:
                if direct_sub_varname is not None:
                    innermost_container[_loc] = param_vars[direct_sub_varname]
                else:
                    innermost_container[_loc] = raw_val.format(**param_vars)
            except KeyError as e:
                msg = f"{e.args[0]} isn't a known parmetrization variable"
                raise ConfReadError(msg, reader, path=path) from None
            except ValueError as e:
                msg = f"issue formatting a parmetrization variable ({e.args[0]})"
                raise ConfReadError(msg, reader, path=path) from None

        if is_scalar:
            obj = obj[0]
        # now, at the very end confirm that we satisfy requirements
        self.req.enforce(obj, reader=reader, path=path)
        return obj


_FORMATTER = string.Formatter()


def _is_template_string(val: TOMLVal, reader: "ConfReader", path: TOMLPath) -> bool:
    any_subs = False
    if isinstance(val, str):
        for _, field_name, _, _ in _FORMATTER.parse(val):
            if field_name is None:
                continue
            field_name = field_name.strip()
            if field_name == "" or field_name[0].isdigit():
                msg = (
                    "unescaped curly brace-pairs must enclose 1 or more "
                    "characters and the first character can't be a digit."
                )
                raise ConfReadError(msg, reader, path=path)
            any_subs = True
    return any_subs


def check_parametrized(
    req: TOMLValReq, reader: "ConfReader", path: TOMLPath
) -> Optional[ParametrizedProxy]:
    path = (path,) if isinstance(path, str) else path

    # holds paths (relative to `path`) corresponding to tables or arrays that we need
    # to visit before we exit this function
    container_stack = []

    def _rel_children_path_list():
        nonlocal container_stack
        rel_path = container_stack.pop()
        _container = reader[path + rel_path]
        if isinstance(_container, list):
            return [rel_path + (i,) for i in range(len(_container))]
        else:
            return [rel_path + (key,) for key in _container]

    first_iter = True
    sub_path_map = {}
    while container_stack or first_iter:
        rel_children_paths = [()] if first_iter else _rel_children_path_list()
        first_iter = False

        for rel_path in rel_children_paths:
            abs_path = path + rel_path
            _tmp = reader[abs_path]

            if _is_template_string(_tmp, reader, path=abs_path):
                sub_path_map[rel_path] = None
            elif isinstance(_tmp, dict) and "use-param" in _tmp:
                _abs_path = path + rel_path + ("use-param",)
                if len(_tmp) > 1:
                    raise ConfReadError("can't have siblings", reader, path=_abs_path)
                sub_path_map[rel_path] = reader[_abs_path]
            elif isinstance(_tmp, (dict, list)):
                container_stack.append(rel_path)
    if len(sub_path_map) != 0:
        return ParametrizedProxy(path, sub_path_map, req)
    return None


class ConfReader:
    # this exists to aggregate related logic and make it easier to provide
    # informative error messages

    file_path: str
    toml_doc: Dict[str, TOMLVal]

    def __init__(self, file_path: os.PathLike):
        self.file_path = str(file_path)
        with open(file_path, "rb") as f:
            self.toml_doc = tomllib.load(f)

    def __repr__(self) -> str:
        return f"ConfReader({self.file_path!r})"

    def __getitem__(self, *keys):
        if len(keys) == 1 and isinstance(keys[0], (tuple, list)):
            keys = keys[0]
        return reduce(operator.getitem, keys, self.toml_doc)

    def get(
        self,
        path: TOMLPath,
        *,
        req: Optional[TOMLValReq] = None,
        allow_parametrized: bool = False,
    ) -> Optional[TOMLVal]:
        """Access the value associated with ``path``

        Parameters
        ----------
        path
            A string, specifying a key in the root table of the document, or a
            sequence of strings where each element is a key of a more nested
            table (analogous to a dotted key in toml).
        req
            Optionally specifies requirements on the parsed value
        allow_parametrize
            Whether the parameter is allowed to be parametrized
        """
        out = self[path]

        tmp = check_parametrized(req, reader=self, path=path)
        if tmp is not None:
            if not allow_parametrized:
                raise ConfReadError("can't be parametrized", self, path)
            return tmp
        elif req is not None:
            req.enforce(val=out, reader=self, path=path)
        return out

    def read_parametrize(
        self, path: TOMLPath
    ) -> List[Dict[str, Union[bool, int, float, str]]]:
        """Attempt to read a parametrize entry from the TOML document

        A parametrize entry is an array of one or more tables, where each
        entry associates keys with one or more values. Currently, the
        keys must be consistent between all tables and values can't be
        arrays or tables.
        """
        # make sure this method functions in a manner consistent with the description
        # provided in PARAMETRIZED_TOML_REQ.description
        #
        # if we need a more sophisticated/powerful mechanism, I think we should mimic
        # the approach of GitHub actions (our current mechanism is equivalent to their
        # strategy.matrix.include parameter). For more info, see
        # https://docs.github.com/en/actions/how-tos/write-workflows/choose-what-workflows-do/run-job-variations

        raw_val = self.get(path)
        if not _is_single_type_toml_arr(raw_val, dict, allow_empty=False):
            raise ConfReadError("isn't an array of 1 or more tables", self, path)

        known_keys = set()
        for i, pack in enumerate(raw_val):
            # check keys
            if len(pack) == 0:
                raise ConfReadError(f"element {i} is an empty table", self, path)
            elif i == 0:
                known_keys = set(pack)
            elif len(known_keys.symmetric_difference(pack)) != 0:
                raise ConfReadError("elements 0 & {i} have different keys", self, path)

            # check values
            for val in pack.values():
                if not isinstance(val, (bool, int, float, str)):
                    msg = (
                        f"element {i} has 1 or more key/value pairs where the value "
                        "isn't a boolean, integer, float, or string"
                    )
                    raise ConfReadError(msg, self, path)
        return raw_val

    def process_table(
        self,
        path: TOMLPath,
        reqs: Collection[TOMLValReq],
    ) -> Tuple[Dict[str, Any], Optional[List[Dict[str, Any]]]]:
        """
        Process a toml table

        Parameters
        ----------
        path
            The path of the table
        reqs
            The requirements to try to parse.

        Returns
        -------
        coerced: dict
            A coerced dict holding values that satisfy all requirements.
        parametrization: list of dicts or None
            If table was parametrized, this is a list of all parametrizations.
        """
        path = (path,) if isinstance(path, str) else path

        # get all keys in the table
        _tab = self[path]
        if not isinstance(_tab, dict):
            raise ConfReadError("isn't a table", self, path=path)
        key_itr = (k for k in _tab.keys() if k != PARAMETRIZED_TOML_REQ.name)
        has_parametrize = PARAMETRIZED_TOML_REQ.name in _tab

        # now let's check for parametrize
        parametrize_path = path + (PARAMETRIZED_TOML_REQ.name,)
        if has_parametrize and (PARAMETRIZED_TOML_REQ.name not in reqs):
            raise ConfReadError("key is forbidden", self, path=parametrize_path)
        elif has_parametrize:
            parametrization = self.read_parametrize(parametrize_path)
        else:
            parametrization = None

        # now, actually parse the table
        coerced = {}
        for key in key_itr:
            full_path = path + (key,)
            req = reqs.get(key)
            if req is None:
                raise ConfReadError("unexpected parameter", self, path=full_path)
            coerced[key] = self.get(
                path=full_path, req=req, allow_parametrized=has_parametrize
            )

        # let's now perform some checks
        if len(reqs) != (len(coerced) + (parametrization is not None)):
            for req in reqs.values():
                if req.name in coerced or not req.required:
                    continue
                msg = "required parameter is missing"
                raise ConfReadError(msg, self, path=path + (req.name,))
        elif parametrization is not None:
            if not any(isinstance(v, ParametrizedProxy) for v in coerced.values()):
                msg = "the parameterization variables are never used"
                raise ConfReadError(msg, self, path=parametrize_path)
        return coerced, parametrization
