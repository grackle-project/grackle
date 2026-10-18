########################################################################
#
# Functions for evolving a fluid container using Grackle
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from collections import defaultdict
import numpy as np

from gracklepy.fluid_container import FluidContainer, _indirectly_calculated_fields

from .physical_constants import \
    gravitational_constant_cgs, \
    sec_per_year

class StateRecorder:
    """
    Aggregate many batches of field-states while evolving a FluidContainer

    The premise is that we can subsequently compute the derived properties
    all at once (for sufficiently complex systems, computing all derived
    properties after a call to fc.solve_chemistry can be more expensive than
    the call to fc.solve_chemistry)
    """

    # number of recordings before object must be reset
    n_recordings_per_reset_: int
    # current number of recordings
    n_recordings_: int

    # number of values copied out of each field when a `FluidContainer` is recorded
    # (if we plan to record a FluidContainer, `fc`, this is `fc.n_vals`)
    field_register_len_: int
    # buffers for recorded field data
    field_data_: dict[str, np.ndarray]
    # tracks recorded times
    time_data_: np.ndarray

    def __init__(
        self,
        n_recordings_per_reset: int,
        field_register_len: int,
        field_names: list[str]
    ):
        # field_names specifies the field names that will be recorded. This should
        # generally just be the input field names (the recorder will NOT compute
        # derived fields)
        assert n_recordings_per_reset > 0
        self.n_recordings_per_reset_ = n_recordings_per_reset
        assert field_register_len > 0
        self.field_register_len_ = field_register_len
        self.n_recordings_ = 0

        field_shape = (n_recordings_per_reset * field_register_len,)
        self.field_data_ = {name: np.empty(shape=field_shape, dtype=np.float64)
                            for name in field_names}
        self.time_data_ = np.empty(shape=field_shape, dtype=np.float64)

    def is_full(self) -> bool:
        return self.n_recordings_ == self.n_recordings_per_reset_

    def is_empty(self) -> bool:
        return self.n_recordings_ == 0

    def record(self, time: float, fc: FluidContainer) -> bool:
        """Record the current state

        This does NOT compute derived fields.
        """
        assert fc.n_vals == self.field_register_len_
        if self.is_full():
            return False

        dst_start = self.n_recordings_ * self.field_register_len_
        dst_stop = dst_start + self.field_register_len_
        for key, dst_buf in self.field_data_.items():
            dst_buf[dst_start:dst_stop] = fc[key]

        self.time_data_[self.n_recordings_] = time

        self.n_recordings_+=1
        return True

    def reset(self):
        self.n_recordings_=0

    def recorded_times_and_fields(self) -> tuple[np.ndarray, dict[str,np.ndarray]]:
        times = self.time_data_[:self.n_recordings_]
        slc = slice(0, self.n_recordings_*self.field_register_len_)
        fields = {k:buf[slc] for k,buf in self.field_data_.items()}
        return (times, fields)

def evolve_freefall(fc, final_density, safety_factor=0.01,
                    include_pressure=True):
    my_chemistry = fc.chemistry_data

    # Set units of gravitational constant
    gravitational_constant = (
        4.0 * np.pi * gravitational_constant_cgs *
        my_chemistry.density_units * my_chemistry.time_units**2)

    # some constants for the analytical free-fall solution
    freefall_time_constant = np.power(((32. * gravitational_constant) /
                                       (3. * np.pi)), 0.5)

    last_n_recorded = 0
    last3_density = np.full(shape=(3, fc.n_vals), fill_value=np.nan, dtype=np.float64)
    last3_pressure = np.full(shape=(3, fc.n_vals), fill_value=np.nan, dtype=np.float64)

    def calc_force_factor(fc):
        nonlocal last_n_recorded
        if include_pressure:
            last3_density[-3:-1, :] = last3_density[-2:, :]
            last3_pressure[-3:-1, :] = last3_pressure[-2:, :]
            fc.calculate_pressure()
            last3_density[-1, :] = fc["density"]
            last3_pressure[-1, :] = fc["pressure"]

            if last_n_recorded < 3:
                last_n_recorded += 1

            fc.calculate_pressure()
            last3_density[-1,:] = fc["density"]
            last3_pressure[-1,:] = fc["pressure"]
            return calculate_collapse_factor(last3_pressure[-last_n_recorded:],
                                             last3_density[-last_n_recorded:])
        else:
            return 0.


    state_recorder = StateRecorder(
        n_recordings_per_reset=100,
        field_register_len=fc.n_vals,
        field_names=fc.input_fields
    )

    data = defaultdict(list)
    current_time = 0.0
    while fc["density"][0] * my_chemistry.density_units < final_density:
        # calculate timestep based on free-fall solution
        dt = safety_factor * \
          np.power(((3. * np.pi) /
                    (32. * gravitational_constant *
                     fc["density"][0])), 0.5)

        state_recorder.record(time=current_time, fc=fc)
        if state_recorder.is_full():
            add_to_data(fc.chemistry_data, state_recorder, data)
            state_recorder.reset() # <- clear out contents

        # compute the new density using the modified
        # free-fall collapse as per Omukai et al. (2005)
        force_factor = calc_force_factor(fc)
        data["force_factor"].append(force_factor)

        # calculate new density from altered free-fall solution
        new_density = np.power((np.power(fc["density"][0], -0.5) -
                                (0.5 * freefall_time_constant * dt *
                                 np.power((1 - force_factor), 0.5))), -2.)

        print("Evolve Freefall - t: %e yr, rho: %e g/cm^3, T: %e K." %
              ((current_time * my_chemistry.time_units / sec_per_year),
               (fc["density"][0] * my_chemistry.density_units),
               fc["temperature"][0]))

        # use this to multiply by elemental densities if you are tracking those
        density_ratio = new_density / fc["density"][0]

        # update densities
        for field in fc.density_fields:
            fc[field] *= density_ratio

        # now update energy for adiabatic heating from collapse
        fc["internal_energy"][0] += (my_chemistry.Gamma - 1.) * \
          fc["internal_energy"][0] * freefall_time_constant * \
          np.power(fc["density"][0], 0.5) * dt

        fc.solve_chemistry(dt)

        # update time
        current_time += dt

    if not state_recorder.is_empty():
        add_to_data(fc.chemistry_data, state_recorder, data)

    for field in data:
        data[field] = np.squeeze(np.array(data[field]))
    return fc.finalize_data(data=data)

def calculate_collapse_factor(pressure, density):
    # Calculate the effective adiabatic index, dlog(p)/dlog(rho).

    if len(pressure) < 3:
        return np.array([0.])

    # compute dlog(p) / dlog(rho) using last two timesteps
    gamma_eff = np.log10(pressure[-1] / pressure[-2]) / \
        np.log10(density[-1] / density[-2])

    # compute a higher order derivative if more than two points available
    if len(pressure) > 2:
        gamma_eff += 0.5 * ((np.log10(pressure[-2] / pressure[-3]) /
                             np.log10(density[-2] / density[-3])) - gamma_eff)

    gamma_eff = np.clip(gamma_eff, a_min=0, a_max=4/3)

    # Equation 9 of Omukai et al. (2005)
    if gamma_eff < 0.83:
        force_factor = gamma_eff * 0
    elif gamma_eff < 1.0:
        force_factor = 0.6 + 2.5 * (gamma_eff - 1) - \
            6.0 * np.power((gamma_eff - 1.0), 2.)
    else:
        force_factor = 1.0 + 0.2 * (gamma_eff - (4./3.)) - \
            2.9 * np.power((gamma_eff - (4./3.)), 2.)
    force_factor = np.clip(force_factor, a_min=0, a_max=0.95)
    return force_factor

def evolve_constant_density(fc, final_temperature=None,
                            final_time=None, safety_factor=0.01):
    my_chemistry = fc.chemistry_data

    if final_temperature is None and final_time is None:
        raise RuntimeError("Must specify either final_temperature " +
                           "or final_time.")

    state_recorder = StateRecorder(
        n_recordings_per_reset=100,
        field_register_len=fc.n_vals,
        field_names=fc.input_fields
    )

    data = defaultdict(list)
    current_time = 0.0
    fc.calculate_cooling_time()
    dt = safety_factor * np.abs(fc["cooling_time"][0])
    fc.calculate_temperature()
    while True:
        if final_temperature is not None and fc["temperature"][0] <= final_temperature:
            break
        if final_time is not None and current_time >= final_time:
            break

        fc.calculate_temperature()
        print("Evolve constant density - t: %e yr, rho: %e g/cm^3, T: %e K." %
              (current_time * my_chemistry.time_units / sec_per_year,
               fc["density"][0] * my_chemistry.density_units,
               fc["temperature"][0]))
        fc.solve_chemistry(dt)

        state_recorder.record(time=current_time, fc=fc)
        if state_recorder.is_full():
            add_to_data(fc.chemistry_data, state_recorder, data)
            state_recorder.reset() # <- clear out contents

        current_time += dt

    if not state_recorder.is_empty():
        add_to_data(fc.chemistry_data, state_recorder, data)

    for field in data:
        data[field] = np.squeeze(np.array(data[field]))
    return fc.finalize_data(data=data)


def add_to_data(chemistry_data, state_recorder, data):
    """
    Add current state_recorder values to the data structure.
    """

    if state_recorder.is_empty():
        return

    times, field_data = state_recorder.recorded_times_and_fields()
    n_vals = field_data["density"].size

    fc = FluidContainer(chemistry_data=chemistry_data, n_vals=n_vals)
    for field_name, buf in field_data.items():
        fc[field_name][:] = buf

    # for indirectly calculated fields, we get more than one field from a call.
    # keep a list of what we've called so we don't have to call them again.
    called = []

    for field in fc.all_fields:
        if field not in fc.input_fields:
            if field in _indirectly_calculated_fields:
                fname = _indirectly_calculated_fields[field]
            else:
                fname = f"calculate_{field}"

            func = getattr(fc, fname, None)
            if func is None:
                raise RuntimeError(f"No function for calculating {field}.")

            if fname not in called:
                func()
                called.append(fname)

        data[field].extend(fc[field])
    data["time"].extend(times)
