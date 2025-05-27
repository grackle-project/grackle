# cimports:
from libcpp cimport bool as CppBool
from libcpp.map cimport map as CppMap
from libcpp.string cimport string as CppString
from libcpp.utility cimport pair
from libcpp.vector cimport vector

from .grackle_defs cimport (
    gr_float,
    c_chemistry_data,
    c_chemistry_data_storage,
    c_code_units,
    c_field_data
)
from .grackle_wrapper cimport GrackleTypePack_ExtType

# runtime imports:
from collections import ChainMap
import functools
import numpy as np
from .fluid_container import FluidContainer, _calculated_fields, _fc_calculated_fields
from .grackle_wrapper import _get_grackle_type_pack
from .utilities.physical_constants import sec_per_year

cdef extern from "grtest/evolve/evolve.hpp" namespace "grtest":
    cdef cppclass IntegrationState:
        IntegrationState()
        unsigned int cycle
        double t

    cdef cppclass EvolveRslt:
        const char* problem_phase
        int task_id
        CppBool triggered_stopping
        EvolveRslt()
        CppString to_err_msg()

    cdef cppclass GrackleTypePack:
        GrackleTypePack()
        c_chemistry_data* my_chem
        c_chemistry_data_storage* my_rates
        c_code_units* my_units
        c_field_data* my_fields

    cdef cppclass IntegratorFn:
        IntegratorFn()
        EvolveRslt operator()(IntegrationState&, unsigned int)

    cdef cppclass IntegratorBuilder:
        IntegratorBuilder()
        IntegratorBuilder& safety_factor(double)
        IntegratorBuilder& max_time(double)
        IntegratorBuilder& min_stopping_temperature_cgs(double)
        IntegratorBuilder& max_stopping_density_cgs(double)
        IntegratorBuilder& config_freefall(double gravitational_constant_cgs,
                                           CppBool pressure_free,
                                           vector[CppString] density_fields,
                                           CppBool strict_density_field_match)
        pair[CppString,IntegratorFn] build(
            const GrackleTypePack& pack,
            CppMap[CppString, gr_float*] chem_register_map,
            CppMap[CppString, gr_float*] output_bufs,
            int buffer_len
        )

_DERIVED_FIELDS = _calculated_fields + _fc_calculated_fields

def _flush_func(cycle, *, flush_cadence, src_dict, dst_dict, preflush_callback):
    """Performs a flush

    Parameters
    ----------
    cycle: int
        The total number of integrated cycles.
    flush_cadence: int
        The nominal number of cycles between flushes
    src_dict: dict
        Data is flushed from this dict
    dst_dict: dict
        Data is flushed to the lists in this dict
    preflush_callback: optional callable
        Called before the flush is performed
    """
    tmp = flush_cadence * (cycle // flush_cadence)
    cycles_since_last_flush = flush_cadence if tmp == cycle else cycle - tmp

    if preflush_callback is not None:
        preflush_callback(cycles_since_last_flush)

    dst_dict["time"].extend(src_dict["time"][0:cycles_since_last_flush])
    register_len = src_dict["density"].size // flush_cadence
    slc = slice(0, register_len * cycles_since_last_flush)
    # we explicitly read the keys from dst_dict (rather than src_dict)
    for key in filter(lambda key: key != "time", dst_dict.keys()):
        dst_dict[key].extend(src_dict[key][slc])


def setup_plumbing(fc, flush_cadence=1, derived_quan=False, extra_buffers=None):
    """Setup plumbing machinery

    Parameters
    ----------
    fc: FluidContainer
        This holds the information used during integration
    flush_cadence: int
        The nominal number of cycles between flushes
    derived_quan: bool
        When True, we computed derived chemistry quantities during flushes.
    extra_buffers
        A sequence of non-chemistry quantity buffers that should be copied.
        This should not contain "time" or any derived quantites

    Returns
    -------
    tmp_buffer: dict
        Used to configure the integrator
    dst_dict: dict
        Where the output data is stored
    flush_fn: callable
        A function used to update dst_dict based on the values in tmp_buffer
    """

    if flush_cadence <= 0:
        raise ValueError("flush_cadence must be positive")
    register_len = fc["density"].size  # number of values updated per cycle
    buffer_len = register_len * flush_cadence
    dtype = fc["density"].dtype

    # most of the work involves assembling tmp_buffers
    tmp_buffers = {"time" : np.empty(buffer_len, dtype=dtype)}

    # add entries for the extra buffers
    if extra_buffers is not None:
        for key in extra_buffers:
            if key == "time" or key in fc or key in _DERIVED_FIELDS:
                raise ValueError(f"extra_buffers can't hold {key}")
            tmp_buffers[key] = np.empty(buffer_len, dtype=dtype)

    # add entries for the input grackle fields. We use FluidContainer as a
    # convenient way to allocate these buffers (it plays a secondary role when
    # derived_quan is True)
    _secondary_fc = FluidContainer(fc.chemistry_data, buffer_len, dtype=dtype)
    for key in _secondary_fc.input_fields:
        tmp_buffers[key] = _secondary_fc[key]

    dst_dict = {key: [] for key in tmp_buffers.keys()}

    # setup _flush_func's kwargs and if relevant, modfiy to dst_dict
    if not derived_quan:
        src_dict_kwarg = tmp_buffers
        callback = None
    else:
        for key in _DERIVED_FIELDS:
            dst_dict[key] = []
        src_dict_kwarg = ChainMap(tmp_buffers, _secondary_fc)

        def callback(cycles_since_last_flush):  # compute derived quantities
            if flush_cadence != cycles_since_last_flush:
                slc = slice(register_len*cycles_since_last_flush, None)
                for key in _secondary_fc.input_fields:
                    _secondary_fc[key][slc] = _secondary_fc[key][0]
            for key in _DERIVED_FIELDS:
                func = getattr(_secondary_fc, f"calculate_{key}")
                assert func is not None, f"No function for calculating {key}"
                func()

    flush_func = functools.partial(
        _flush_func,
        flush_cadence=flush_cadence,
        src_dict=src_dict_kwarg,
        dst_dict=dst_dict,
        preflush_callback=callback
    )

    return tmp_buffers, dst_dict, flush_func


cdef CppMap[CppString, gr_float*] construct_ptr_map(object d, object keys=None):
    cdef CppMap[CppString, gr_float*] out
    cdef CppString tmp_string
    cdef gr_float[::1] tmp_view
    cdef object it = d.keys() if keys is None else keys
    for key in it:
        tmp_string = key.encode("ascii")
        tmp_view = d[key]
        out[tmp_string] = <gr_float *> &tmp_view[0]
    return out


cdef object invoke_integrator(
    object fc,
    const IntegratorBuilder& builder,
    unsigned int flush_cadence=1,
    object derived_quan=False,
    object extra_buffers=None,
    object log_prefix=""
) except+:
    # Constructs an integrator, use it to evolve the input data, and return a
    # dictionary of arrays holding the values computed during each integration
    # cycle
    #
    # Parameters
    # ----------
    # fc: The fluid container to evolve
    # builder: contains configuration information about the integrator
    # flush_cadence: The nominal number of cycles between flushes
    # derived_quan: whether to compute derived chemistry quantities
    # extra_buffers: A sequence of non-chemistry quantity buffers that should
    #    be copied into the output array. This should not contain "time" or any
    #    derived quantites. This will contain quantities like "force_factor"
    # log_prefix: A prefix for each log message

    # Step 1: Prepare to build the integrator
    cdef GrackleTypePack_ExtType pack_ext_type = _get_grackle_type_pack(fc)
    cdef GrackleTypePack pack
    pack.my_chem = pack_ext_type.my_chem
    pack.my_rates = pack_ext_type.my_rates
    pack.my_units = pack_ext_type.my_units
    pack.my_fields = &pack_ext_type.my_fields

    tmp_buffers, data, flush_func = setup_plumbing(
        fc,
        flush_cadence=flush_cadence,
        derived_quan=derived_quan,
        extra_buffers=extra_buffers
    )

    # Step 2: construct the integrator
    cdef pair[CppString, IntegratorFn] tmp_pair = builder.build(
        pack=pack,
        # in the future, chem_register_map won't be necessary
        chem_register_map=construct_ptr_map(fc, keys=fc.input_fields),
        output_bufs=construct_ptr_map(tmp_buffers),
        buffer_len=tmp_buffers["density"].size
    )

    if (tmp_pair.first.size() != 0):
        raise RuntimeError(
            "Error while creating integrator: " + tmp_pair.first.decode("ascii")
        )
    cdef IntegratorFn integrator = tmp_pair.second

    # Step 3: invoke the integrator
    cdef IntegrationState state  # default constructor sets cycle=0 & t =0
    cdef EvolveRslt tmp_rslt  # forward declaration
    cdef object time_units = fc.chemistry_data.time_units
    cdef object density_units = fc.chemistry_data.density_units

    while True:
        tmp_rslt = integrator(state, flush_cadence)
        if tmp_rslt.problem_phase != NULL:
            raise RuntimeError(tmp_rslt.to_err_msg().decode("ascii"))
        flush_func(state.cycle) # flush the data
        t = data["time"][-1] * time_units / sec_per_year
        rho = data["density"][-1] * density_units
        T = data["temperature"][-1]
        print(f"{log_prefix}t: {t:e} yr, rho: {rho:e} g/cm^3, T: {T:e} K.")
        if tmp_rslt.triggered_stopping:
            return data

def evolve_constant_density(
    fc, final_temperature = None, final_time = None, safety_factor = 0.01,
):
    """Evolves the fluid_container while holding density constant
    """
    cdef IntegratorBuilder builder

    if final_temperature is not None:
        # the value is checked while building the integrator
        builder.min_stopping_temperature_cgs(final_temperature)

    if final_time is not None:
        builder.max_time(final_time)

    builder.safety_factor(safety_factor)

    # determine flux_cadence (this is fairly ad hoc). Considerations:
    # - we bound flux_cadence based on register_len since we need to allocate
    #   register_len * flux_cadence * len(fc.input_fields) floating point vals
    # - from some quick and dirty test, I saw marginal benefits for values
    #   above 20 and register_len == 1
    # - maybe we should make this a kwarg in the future?
    register_len = fc["density"].size
    flush_cadence = max(20, min(256 // register_len, 1))

    data = invoke_integrator(
        fc=fc,
        builder=builder,
        flush_cadence=flush_cadence,
        derived_quan=True,  # <- for historical purposes
        extra_buffers=[],
        log_prefix="Evolve constant density - "
    )
    for field in data:
        data[field] = np.array(data[field])
    return fc.finalize_data(data=data)
