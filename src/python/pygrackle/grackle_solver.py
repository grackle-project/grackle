########################################################################
#
# Implements the GrackleSolver class
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from .fluid_container import FluidContainer, _expected_density_fields
from .grackle_wrapper import (
    chemistry_data,
    _wrapped_c_chemistry_data,
    solve_chemistry,
    calculate_property
)
from .utilities.arraymap import ArrayMap



_ALWAYS_NEEDED_UNIT_ATTRS = ("density_units", "length_units", "time_units")
_UNITS_ATTRS = frozenset(_ALWAYS_NEEDED_UNIT_ATTRS +
                         ("velocity_units", "a_value", "a_units",
                          "comoving_coordinates"))

_CHEM_PARAMETERS = frozenset(_wrapped_c_chemistry_data.int_keys() +
                             _wrapped_c_chemistry_data.double_keys() +
                             _wrapped_c_chemistry_data.string_keys())

def _is_rateattr(attr, descr):
    # function that is passed the name of a chemistry_data attribute and an
    # object that may implement the "descriptor protocol" (e.g. a property or
    # method)
    # -> these arguments are the key-value pairs in chemistry_data.__dict__
    # -> this function returns True if the arguments specifically refer to a
    #    rates data attribute

    if len(attr) == 0 or attr[0] == '_':
        return False # ignore private attributes/magic methods
    elif (attr in _UNITS_ATTRS) or (attr in _CHEM_PARAMETERS):
        # ignore attributes associated with code_units struct & _wrapped_c_chemistry_data
        return False

    # https://docs.python.org/3/howto/descriptor.html#descriptor-protocol says
    # that all descriptors define __get__, __set__ or __delete__ methods. It
    # also says that "data descriptors" (like properties) always define a
    # __set__ or __delete__ method (read-only data descriptors define __set__,
    # but have it raise AttributeError)
    is_datadescr = hasattr(descr, '__set__') or hasattr(descr, '__delete__')
    return hasattr(descr, '__get__') and is_datadescr


_RATES_ATTRS = frozenset(
    k for (k,v) in chemistry_data.__dict__.items() if _is_rateattr(k, v))



# All of this machinery here is designed to help us provide a nice interface
# for querying/mutating chemistry-parameters and rates:
#        `print(grackle_solver.parameters.primordial_chemistry)`
#        `print(grackle_solver.rates.k13)`
#        `print(grackle_solver.code_units.a_value)`
# if we decide we like this, we can clean things up inside the cython file to
# implement this all more directly.
# -> For example, `grackle_solver.parameters` might directly return an instance
#    of the class that we now call `_wrapped_c_chemistry_data` (we would need
#    to refactor that class to make it look more like a namedtuple/dataclass
#    than a dictionary, but I think that makes sense since there are typing
#    requirements)
# -> we could introduce extension types similar to `_wrapped_c_chemistry_data`
#    that wrap the `code_units` and `chemistry_data_storage` structs
# Other questions:
# -> how do we feel about mutability? I like the idea that
#    `grackle_solver.code_units` provides an immutable object (we already have
#    the user specify `a_value` as an argument to the relevant functions)
# -> I also kinda like the idea of making `grackle_solver.parameters` produce
#    an immutable object, but in the short-term, I can see reasons not to do
#    that...


class _DataProxy:
    """
    A proxy class for accessing (a subset of) attributes of a wrapped class.

    It would be nice to provide a dataclass-like interface
    """

    # we intentionally use __slots__ to help users avoid typos when they are
    # trying to modify a value
    __slots__ = ("_wrapped_attrs", "_wrapped_obj",  "_mutable")

    def __init__(self, wrapped_attrs, wrapped_obj, mutable = True):
        self._wrapped_attrs = wrapped_attrs # do this first!
        self._wrapped_obj = wrapped_obj
        self._mutable = mutable

    def __getattribute__(self, name):
        wrapped_attrs = object.__getattribute__(self, "_wrapped_attrs")
        if name in wrapped_attrs:
            wrapped_obj = object.__getattribute__(self, "_wrapped_obj")
            return getattr(wrapped_obj, name)
        else:
            return object.__getattribute__(self, name)

    def __setattr__(self, name, value):
        if name == "_wrapped_attrs":
            return object.__setattr__(self, name, value)

        wrapped_attrs = object.__getattribute__(self, "_wrapped_attrs")
        if name in wrapped_attrs:
            if object.__getattribute__(self, "_mutable"):
                wrapped_obj = object.__getattribute__(self, "_wrapped_obj")
                return setattr(wrapped_obj, name, value)
            raise AttributeError("Immutable proxies can't mutate the wrapped "
                                 "object's attributes")
        else:
            return object.__setattr__(self, name, value)

    def asdict(self):
        """
        Returns a dictionary of name,value pairs of the attributes accessible
        in the underlying class

        Note
        ----
        The name is based on the function name used with dataclasses
        """
        return dict((attr, getattr(self, attr)) for attr in self._wrapped_attrs)


def _configure_codeunits(units_obj, code_units_kwargs):
    """
    Setup values related to the code_units struct.

    Notes
    -----
    It would be easier to make an extension type that wraps code_units and make
    the constructor responsible for this sort of error checking
    """

    units_obj.comoving_coordinates = code_units_kwargs.get(
        "comoving_coordinates", 0)
    if units_obj.comoving_coordinates not in [0, 1]:
        raise ValueError("'comoving_coordinates' code_units kwarg must be "
                         "associated with a value of 0 or 1")

    if "a_value" in code_units_kwargs:
        units_obj.a_value = code_units_kwargs["a_value"]
        if "a_units" not in code_units_kwargs:
            raise ValueError("'a_value' code_units kwarg cannot be specified "
                             "without the 'a_units' kwarg")
        units_obj.a_units = code_units_kwargs["a_units"]
    elif "a_units" in code_units_kwargs:
        raise ValueError("'a_units' code_units kwargs cannot be specified "
                         "without the 'a_value' kwarg")
    elif units_obj.comoving_coordinates == 0:
        units_obj.a_value = 0
        units_obj.a_units = 0
    else:
        raise ValueError("when the 'comoving_coordinates' code_units kwarg "
                         "has a value of 1, the 'a_value' and 'a_units' "
                         "code_units key-value pairs must be specified")

    for kw in _ALWAYS_NEEDED_UNIT_ATTRS:
        if kw not in code_units_kwargs:
            raise ValueError(f"'{kw}' code_units kwarg must be specified")
        setattr(units_obj, kw, code_units_kwargs[kw])

    vel_u = code_units_kwargs.get("velocity_units", None)
    if vel_u is None:
        units_obj.set_velocity_units()
    else:
        units_obj.velocity_units = vel_u

    # last check:
    for key in code_units_kwargs:
        if key not in _UNITS_ATTRS:
            raise ValueError(f"'{key}' is not a known code_units kwarg")

def _get_fc(grackle_solver, array_map = None):
    if array_map is None:
        return FluidContainer(grackle_solver._pychemistry_data, n_vals = 1)
    else:
        assert isinstance(array_map, ArrayMap)
        array_shape = array_map.array_shape
        assert len(array_shape) == 1
        # direclty overwrite the arrays in newly allocaed FluidContainer with
        # arrays allocated in array_map
        fc = FluidContainer(grackle_solver._pychemistry_data,
                            n_vals = array_shape[0])
        # FluidContainer holds keys for derived quantities in addition to
        # required quantities
        for k in array_map:
            fc[k] = array_map[k]
        return fc

class GrackleSolver:
    """
    Represents a fully initialized grackle-chemistry solver

    Parameters
    ----------
    chemistry_parameters : dict
        Holds kwargs that directly map to the members of the chemistry_data
        struct provided by grackle
    code_units : dict
        Holds kwargs that directly map to the members of the chemistry_data
        struct provided by grackle
    """

    def __init__(self, chemistry_parameters, code_units):
        self._pychemistry_data = chemistry_data._from_wrapped_c_chemistry_data(
            _wrapped_c_chemistry_data(chemistry_parameters)
        )
        _configure_codeunits(self._pychemistry_data, code_units)
        tmp = self._pychemistry_data.initialize()
        if tmp != 1:
            raise ValueError("there was an issue initializing grackle-solver")
        # from here on out, we NEVER mutate self._pychemistry_data again!

    def expected_species_density_fields(self):
        """
        Return the required number-density fields
        """
        chem_data = self._pychemistry_data
        tmp = _expected_density_fields(
            primordial_chemistry = chem_data.primordial_chemistry,
            metal_cooling = chem_data.metal_cooling,
            use_dust_density_field = chem_data.use_dust_density_field)
        if 'density' in tmp:
            tmp.remove('density')
        return tmp

    @property
    def parameters(self):
        """
        Returns read-only proxy that provides access to parameters information
        that was used to configure the `GrackleSolver`

        If you have a `GrackleSolver` instance called `gs`, you might use
        `print(gs.parameters.primordial_chemistry)` to print the value of
        `primordial_chemistry`. You might call `gs.code_units.asdict()` to get
        a dictionary containing all relevant accessible information.

        Notes
        -----
        The user should generally treat this as immutable (even if we
        temporarily allow mutation). Modifying values after initialization can
        cause Grackle to find itself in an unexpected state.
        """
        return _DataProxy(_CHEM_PARAMETERS, self._pychemistry_data,
                          mutable = True)

    @property
    def code_units(self):
        """
        Returns read-only proxy that provides access to units-information that
        was used to configure the `GrackleSolver`

        If you have a `GrackleSolver` instance called `gs`, you might use
        `print(gs.code_units.a_value)` to print the value of `a_value`. You
        might call `gs.code_units.asdict()` to get a dictionary containing all
        relevant accessible information.

        Notes
        -----
        Under `GrackleSolver`'s current design, this MUST return an immutable
        object because the current design does not provide the user with any
        way to automatically reset the velocity units after modifying other
        related quantities.

        In reality, the user shouldn't have any need to change any unit value
        after initializing a `GrackleSolver` instance other than `a_value`, and
        `a_value` can be temporarily modified in all relevant methods by using
        the `a_value` kwarg.
        """
        return _DataProxy(_UNITS_ATTRS, self._pychemistry_data, mutable = False)

    @property
    def rates(self):
        """
        Returns proxy that provides access to rates information that was
        initialized during the configuration of the `GrackleSolver`

        If you have a `GrackleSolver` instance called `gs`, you might use
        `print(gs.rates.k6)` to print the value of `gs.rates.k6`. You can also
        update the value of `gs.rates.k6`. You might call `gs.rates.asdict()`
        to get a dictionary containing all relevant accessible information.
        """
        return _DataProxy(_RATES_ATTRS, self._pychemistry_data, mutable = True)

    def solve_chemistry(self, data, dt, *, current_a_value = None,
                        inplace = False, grid_dx = None):
        """
        Evolves the species densities and internal energies over a given
        timestep by solving the chemistry and cooling rate equations

        Parameters
        ----------
        data : mapping
            Mapping (e.g. a `dict`) that associates field names (required by
            grackle) with ordinary 1D numpy arrays (i.e. they can't currently
            be subclasses), that share a common length. The values of the field
            are understood to be in code-units.
        dt : float
            The integration timestep (in code-units)
        current_a_value : float, optional
            Optional parameter to overwrite the cosmological scale factor used
            in the calculation.
        inplace : bool, optional
            When `True`, this function attempts to update the field values
            in-place. By default, a copy of the fields are mutated.
        grid_dx : float, optional
            Optionally specifies the cell-width. This is only used for certain
            kinds of self-shielding methods.

        Returns
        -------
        mapping
            Mapping containing the fields after they have been evolved. When
            `inplace` is `True`, the mapping may not match the one used
            specified in `data`, the contained arrays will directly reference
            the arrays in `data`. This will also contain fields used in the
            calculation that are not necessarily updated.
        """
        # currently, grid_dx is just a placeholder. Eventually, we will pass it
        # through solve_chemistry (once we add support for 3D arrays)
        if grid_dx is not None:
            raise NotImplementedError(
                "We currently don't support a non-None grid_dx value")

        array_map = ArrayMap(data, copy = inplace)
        fc = _get_fc(self, array_map = array_map)
        solve_chemistry(fc, dt, current_a_value = current_a_value)
        return array_map

    def calculate_cooling_time(self, data, *, current_a_value = None,
                               inplace = False, grid_dx = None):
        """
        Calculate the instantaneous cooling time.

        Parameters
        ----------
        data : mapping
            Mapping (e.g. a `dict`) that associates field names (required by
            grackle) with ordinary 1D numpy arrays (i.e. they can't currently
            be subclasses), that share a common length. The values of the field
            are understood to be in code-units.
        current_a_value : float, optional
            Optional parameter to overwrite the cosmological scale factor used
            in the calculation.
        inplace : bool, optional
            By default, a copy of the fields is used. When `True`, this
            function performs the calculation by directly using the arrays held
            by `data`. To reduce memory overhead, calculation may involve
            modifying field values (e.g. for unit-conversions). While the
            function includes logic to the field values back to their original
            values, the final state of the arrays may not exactly match the
            original state (due to the nature of floating point operations and,
            in some cases, the application of floors).
        grid_dx : float, optional
            Optionally specifies the cell-width. This is only used for certain
            kinds of self-shielding methods.

        Returns
        -------
        np.ndarray
            Array with the same shape as the arrays in `data` that holds the
            instantaneous cooling_time (in code units). Negative values
            correspond to cooling and positive values correspond to heating.
        """
        # currently, grid_dx is just a placeholder. Eventually, we will pass it
        # through calculate_property (once we add support for 3D arrays)
        if grid_dx is not None:
            raise NotImplementedError(
                "We currently don't support a non-None grid_dx value")

        return calculate_property(
            fc = _get_fc(self, array_map = ArrayMap(data, copy = inplace)),
            property_name = "cooling_time", current_a_value = current_a_value)

    def calculate_dust_temperature(self, data, *, current_a_value = None,
                                   inplace = False):
        """
        Calculates the dust temperature.

        Parameters
        ----------
        data : mapping
            Mapping (e.g. a `dict`) that associates field names (required by
            grackle) with ordinary 1D numpy arrays (i.e. they can't currently
            be subclasses), that share a common length. The values of the field
            are understood to be in code-units.
        current_a_value : float, optional
            Optional parameter to overwrite the cosmological scale factor used
            in the calculation.
        inplace : bool, optional
            By default, a copy of the fields is used. When `True`, this
            function performs the calculation by directly using the arrays held
            by `data`. To reduce memory overhead, calculation may involve
            modifying field values (e.g. for unit-conversions). While the
            function includes logic to the field values back to their original
            values, the final state of the arrays may not exactly match the
            original state (due to the nature of floating point operations and,
            in some cases, the application of floors).

        Returns
        -------
        np.ndarray
            Array with the same shape as the arrays in `data` that holds the
            dust temperature (in units of Kelvin).
        """
        return calculate_property(
            fc = _get_fc(self, array_map = ArrayMap(data, copy = inplace)),
            property_name = "dust_temperature",
            current_a_value = current_a_value)

    def calculate_gamma(self, data, *, current_a_value = None,
                        inplace = False):
        """
        Calculate the adiabatic index. This is only useful in configurations
        that track molecular Hydrogen.

        Parameters
        ----------
        data : mapping
            Mapping (e.g. a `dict`) that associates field names (required by
            grackle) with ordinary 1D numpy arrays (i.e. they can't currently
            be subclasses), that share a common length. The values of the field
            are understood to be in code-units.
        current_a_value : float, optional
            Optional parameter to overwrite the cosmological scale factor used
            in the calculation.
        inplace : bool, optional
            By default, a copy of the fields is used. When `True`, this
            function performs the calculation by directly using the arrays held
            by `data`. To reduce memory overhead, calculation may involve
            modifying field values (e.g. for unit-conversions). While the
            function includes logic to the field values back to their original
            values, the final state of the arrays may not exactly match the
            original state (due to the nature of floating point operations and,
            in some cases, the application of floors).

        Returns
        -------
        np.ndarray
            Array with the same shape as the arrays in `data` that holds the
            adiabatic index.
        """
        return calculate_property(
            fc = _get_fc(self, array_map = ArrayMap(data, copy = inplace)),
            property_name = "gamma", current_a_value = current_a_value)

    def calculate_pressure(self, data, *, current_a_value = None,
                           inplace = False):
        """
        Calculates the thermal gas pressure.

        Parameters
        ----------
        data : mapping
            Mapping (e.g. a `dict`) that associates field names (required by
            grackle) with ordinary 1D numpy arrays (i.e. they can't currently
            be subclasses), that share a common length. The values of the field
            are understood to be in code-units.
        current_a_value : float, optional
            Optional parameter to overwrite the cosmological scale factor used
            in the calculation.
        inplace : bool, optional
            By default, a copy of the fields is used. When `True`, this
            function performs the calculation by directly using the arrays held
            by `data`. To reduce memory overhead, calculation may involve
            modifying field values (e.g. for unit-conversions). While the
            function includes logic to the field values back to their original
            values, the final state of the arrays may not exactly match the
            original state (due to the nature of floating point operations and,
            in some cases, the application of floors).

        Returns
        -------
        np.ndarray
            Array with the same shape as the arrays in `data` that holds the
            thermal gas pressure (in code units).
        """
        return calculate_property(
            fc = _get_fc(self, array_map = ArrayMap(data, copy = inplace)),
            property_name = "pressure", current_a_value = current_a_value)

    def calculate_temperature(self, data, *, current_a_value = None,
                              inplace = False):
        """
        Calculates the gas temperature.

        Parameters
        ----------
        data : mapping
            Mapping (e.g. a `dict`) that associates field names (required by
            grackle) with ordinary 1D numpy arrays (i.e. they can't currently
            be subclasses), that share a common length. The values of the field
            are understood to be in code-units.
        current_a_value : float, optional
            Optional parameter to overwrite the cosmological scale factor used
            in the calculation.
        inplace : bool, optional
            By default, a copy of the fields is used. When `True`, this
            function performs the calculation by directly using the arrays held
            by `data`. To reduce memory overhead, calculation may involve
            modifying field values (e.g. for unit-conversions). While the
            function includes logic to the field values back to their original
            values, the final state of the arrays may not exactly match the
            original state (due to the nature of floating point operations and,
            in some cases, the application of floors).

        Returns
        -------
        np.ndarray
            Array with the same shape as the arrays in `data` that holds the
            gas temperature (in units of Kelvin).
        """
        return calculate_property(
            fc = _get_fc(self, array_map = ArrayMap(data, copy = inplace)),
            property_name = "temperature", current_a_value = current_a_value)
