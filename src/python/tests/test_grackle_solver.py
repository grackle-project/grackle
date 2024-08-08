########################################################################
#
# Tests the GrackleSolver class. To start, most of the tests compare
# results with the older API
#
########################################################################

import os
import numpy as np
import pytest

from pygrackle import FluidContainer, chemistry_data
from pygrackle.fluid_container import _expected_density_fields
from pygrackle.grackle_solver import GrackleSolver
from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc, \
    boltzmann_constant_cgs
from pygrackle.utilities.testing import assert_allclose_resultdict
from pygrackle.utilities.units import set_cosmology_units

def _fields_updated_by_chemistry_solve(chemistry_obj):
    # first get list of relevant density-related fields
    if isinstance(chemistry_obj, GrackleSolver):
        tmp = chemistry_obj.expected_species_density_fields()
    elif isinstance(chemistry_obj, chemistry_data):
        # we avoid using FluidContainer's density_fields property since the
        # property may include unnecessary fields (e.g. the "dust" field in
        # cases where the "dust" density field is unused)
        tmp = _expected_density_fields(
            primordial_chemistry = chemistry_obj.primordial_chemistry,
            metal_cooling = chemistry_obj.metal_cooling,
            use_dust_density_field = chemistry_obj.use_dust_density_field)
    else:
        raise TypeError("chemistry_obj has unexpected type")
    # now modify that list and return result
    return [e for e in tmp if e not in {"density", "metal"}] + ["energy"]


_derived_prop_pairs = [
    ('cooling_time', 'calculate_cooling_time'),
    ('gamma', 'calculate_gamma'),
    ('pressure', 'calculate_pressure'),
    ('temperature', 'calculate_temperature'),
    ('dust_temperature', 'calculate_dust_temperature')]

def _perform_calculations(chemistry_obj, initial_fluid_config,
                          nominal_dt = None, current_a_value = None,
                          inplace = False):
    """
    Perform standard battery of tests
    """
    results = {} # dictionary containing the results

    if nominal_dt is None:
        # fall back to a timestep of 1 kyr in proper units (no idea if this is
        # actually reasonable or not)
        if isinstance(chemistry_obj, GrackleSolver):
            _time_units = chemistry_obj.code_units.time_units
        else:
            _time_units = chemistry_obj.time_units
        dt = 1000 * 3.156e7 / _time_units
    else:
        dt = nominal_dt

    # first, initialize objects that hold the input data & define the functions
    # that reset contents with initial values
    if isinstance(chemistry_obj, GrackleSolver):
        # when inplace == False, we should not need to ever reset the values
        # after the first initialization of data_dict
        reset_chem_state = lambda: None
        data_dict = {}

        def reset_initial_field_data():
            for (k,arr) in initial_fluid_config.items():
                data_dict[k] = arr.copy()
    else:
        orig_a_value = chemistry_obj.a_value

        def reset_chem_state():
            chemistry_obj.a_value = orig_a_value

        if current_a_value is None:
            current_a_value = orig_a_value
        fc = FluidContainer(chemistry_obj,
                            n_vals = initial_fluid_config['density'].size)

        def reset_initial_field_data():
            # if we try overwrite the contents of fields in fc that don't have
            # corresponding entries in initial_fluid_config so that the fields
            # are uniformly 0, weird issues arise
            for k in initial_fluid_config:
                fc[k][:] = initial_fluid_config[k]

    # next, compute derived properties
    for field, fn_name in _derived_prop_pairs:
        reset_initial_field_data()  # <- reset initial state
        if isinstance(chemistry_obj, GrackleSolver):
            rslt = getattr(chemistry_obj,fn_name)(
                data_dict, current_a_value = current_a_value, inplace = inplace)
        else:
            fc.chemistry_data.a_value = current_a_value
            rslt = fc[field]
        results[f'derived-{fn_name}/{field}'] = rslt

    # now, compute the values after executing the chemistry-solver
    reset_initial_field_data()  # <- reset initial state
    if isinstance(chemistry_obj, GrackleSolver):
        chem_solve_rslt = chemistry_obj.solve_chemistry(
            data_dict, dt = dt, current_a_value = current_a_value,
            inplace = inplace)
    else:
        fc.chemistry_data.a_value = current_a_value
        fc.solve_chemistry(dt = dt)
        chem_solve_rslt = fc
    # copy the results from solve_chemistry
    for field in _fields_updated_by_chemistry_solve(chemistry_obj):
        results[f'post-solve_chemistry/{field}'] = chem_solve_rslt[field]

    reset_chem_state()
    return results

def _simple_initial_vals(chemistry_parameters, code_units, *,
                         initial_density_cgs = 2e-25,
                         nominal_initial_temperature_cgs = 1.e6,
                         nominal_initial_mmw = 0.6,
                         nominal_initial_gamma = 5/3.0,
                         nominal_HMassFrac = 0.76,
                         nominal_SolarMetalFractionByMass = 0.0125):
    # We intentionally try to avoid using internal grackle functionality in
    # this particular function. We use the term, nominal, to reflect the fact
    # that the values may not necessarily reflect the actual values that are
    # being used

    tiny_number = 1e-20

    make_array = lambda val: np.full((1,), dtype = 'f8', fill_value = val)
    out = {}
    out["density"] = make_array(initial_density_cgs /
                                code_units['density_units'])

    if chemistry_parameters['primordial_chemistry'] > 0:
        out["HI"] = nominal_HMassFrac * out["density"]
        out["HII"] = tiny_number * out["density"]
        out["HeI"] = (1.0 - nominal_HMassFrac) * out["density"]
        out["HeII"] = tiny_number * out["density"]
        out["HeIII"] = tiny_number * out["density"]
        out["de"] = tiny_number * out["density"]
    if chemistry_parameters['primordial_chemistry'] > 1:
        out["H2I"] = tiny_number * out["density"]
        out["H2II"] = tiny_number * out["density"]
        out["HM"] = tiny_number * out["density"]
    if chemistry_parameters['primordial_chemistry'] > 2:
        out["DI"] = 2.0 * 3.4e-5 * out["density"]
        out["DII"] = tiny_number * out["density"]
        out["HDI"] = tiny_number * out["density"]
    if chemistry_parameters['metal_cooling'] == 1:
        out["metal"] = 0.1 * out["density"] * nominal_SolarMetalFractionByMass

    out["x-velocity"] = 0.0 * out["density"]
    out["y-velocity"] = 0.0 * out["density"]
    out["z-velocity"] = 0.0 * out["density"]

    nominal_initial_ndens_cgs = initial_density_cgs / (mass_hydrogen_cgs *
                                                       nominal_initial_mmw)
    nominal_initial_pressure_cgs = (
        nominal_initial_ndens_cgs * boltzmann_constant_cgs *
        nominal_initial_temperature_cgs)

    initial_energy_cgs = (
        nominal_initial_pressure_cgs / ((nominal_initial_gamma - 1.0) *
                                        initial_density_cgs))
    energy_units = (code_units["length_units"] / code_units["time_units"]) ** 2

    out["energy"] = make_array(initial_energy_cgs / energy_units)
    return out

def _initial_props(primordial_chemistry = 0, initial_redshift = 0.5,
                   comoving = False):

    my_dir = os.path.dirname(os.path.abspath(__file__))

    chemistry_parameters = dict(
        use_grackle = 1,
        with_radiative_cooling = 1,
        primordial_chemistry = primordial_chemistry,
        metal_cooling = 1,
        UVbackground = 1,
        self_shielding_method = 0,
        H2_self_shielding = 0,
        grackle_data_file = os.path.join(
            my_dir, "..", "..", "..", "input", "CloudyData_UVB=HM2012.h5")
    )

    if comoving:
        code_units = set_cosmology_units(None,
                                         current_redshift = initial_redshift,
                                         initial_redshift = initial_redshift)
    else:
        a_units = 1.0
        code_units = dict(
            comoving_coordinates = 0, # proper units
            a_units = a_units,
            a_value = 1. / (1. + initial_redshift) / a_units,
            density_units = mass_hydrogen_cgs, # rho = 1.0 is 1.67e-24 g
            length_units = cm_per_mpc,         # 1 Mpc in cm
            time_units = sec_per_Myr,          # 1 Myr in s
        )

    initial_fluid_config = _simple_initial_vals(chemistry_parameters,
                                                code_units)
    nominal_dt = 1000 * 3.156e7 / code_units['time_units'] # <- 1 kyr

    return chemistry_parameters, code_units, initial_fluid_config, nominal_dt

@pytest.mark.parametrize("current_a_value", [None, 0.75])
@pytest.mark.parametrize("primordial_chemistry", [0,1,2,3])
def test_API_comparison(primordial_chemistry, current_a_value):
    initial_redshift, comoving = 0.5, False
    if current_a_value is not None:
        comoving = True
        initial_a_value = (1.0/ (1.0 + initial_redshift))
        if initial_a_value > current_a_value:
            raise ValueError(
                "Invalid choice for current_a_value. It should be None or "
                f"it must correspond to a redshift after {initial_redshift}")

    tmp = _initial_props(primordial_chemistry = primordial_chemistry,
                         initial_redshift = initial_redshift,
                         comoving = comoving)
    chemistry_parameters, code_units, initial_fluid_config, nominal_dt = tmp

    rslt_pack = {}
    for approach in ('FluidContainer', 'GrackleSolver',
                     'GrackleSolver-inplace'):
        if approach == 'FluidContainer':
            my_chemistry = chemistry_data()
            for k,v in chemistry_parameters.items():
                setattr(my_chemistry, k, v)
            for k,v in code_units.items():
                setattr(my_chemistry, k, v)
            my_chemistry.set_velocity_units()
            if my_chemistry.initialize() != 1:
                raise RuntimeError("There was an issue with initialization")
            chemistry_obj = my_chemistry
            inplace = None
        elif 'GrackleSolver' == approach:
            chemistry_obj = GrackleSolver(chemistry_parameters, code_units)
            inplace = False
        elif 'GrackleSolver-inplace' == approach:
            chemistry_obj = GrackleSolver(chemistry_parameters, code_units)
            inplace = True
        else:
            raise RuntimeError()

        rslt_pack[approach] = _perform_calculations(
            chemistry_obj, initial_fluid_config,
            current_a_value = current_a_value,
            nominal_dt = nominal_dt, inplace = inplace)

    for (key, descr) in [('GrackleSolver', 'are not'),
                         ('GrackleSolver-inplace', 'are')]:
        assert_allclose_resultdict(
            rslt_pack[key], reference = rslt_pack['FluidContainer'],
            cmp_kwarg_spec = dict(), # require exact matches
            err_msg=('The reference values use the FluidContainer API while '
                     'the actual values used the GrackleSolver API (the '
                     f'arrays {descr} directly used inplace)')
        )

