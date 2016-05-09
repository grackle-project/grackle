########################################################################
#
# Units helper functions
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np

def set_cosmology_units(my_units, hubble_constant=0.704,
                        omega_matter=0.268, omega_lambda=0.732,
                        current_redshift=0.0, initial_redshift=0.0,
                        comoving_box_size=1.0):
    """
    Set cosmological units like Enzo.

    Greg Bryan's note on cosmology units:

    time:        utim = 1 / sqrt(4 * \pi * G * \rho_0 * (1+zri)^3)
    density:     urho = \rho_0 * (1+z)^3
    length:      uxyz = (1 Mpc) * box / h / (1+z)
    velocity:    uvel = uaye * (uxyz / a) / utim  (since u = a * dx/dt)
    (*)  temperature: utem = m_H * \mu / k * uvel**2
    a(t):        uaye = 1 / (1 + zri)

    where:
    box     - size of simulation box in Mpc/h
    zri     - initial redshift (start of simulation)
    \rho_0  = 3*\Omega_0*H_0^2/(8*\pi*G)
    Omega_0 - the fraction of non-relativistic matter at z=0

    Note that two definitions are dependent on redshift (urho
    and uxyz) so make sure to call this routine immediately
    before writing.

    * - the utem given below assumes that \mu = 1, so you must
    multiply the resulting temperature field by \mu.
    """

    my_units.comoving_coordinates = 1
    my_units.a_units = 1. / (1. + initial_redshift)
    my_units.a_value = 1.0 / (1.0 + current_redshift) / \
        my_units.a_units
    my_units.density_units = 1.8788e-29 * omega_matter * \
        np.power(hubble_constant, 2) * np.power(1 + current_redshift, 3)
    my_units.length_units = 3.085678e24 * comoving_box_size / \
        hubble_constant / (1. + current_redshift)
    my_units.time_units = 2.519445e17 / np.sqrt(omega_matter) / \
        hubble_constant / np.power(1 + initial_redshift, 1.5)
    my_units.velocity_units = 1.22475e7 * comoving_box_size * \
        np.sqrt(omega_matter) * np.sqrt(1 + initial_redshift)
