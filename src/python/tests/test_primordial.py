########################################################################
#
# Chemistry testing functions
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

from matplotlib import pyplot

from pygrackle.utilities.testing import \
    assert_rel_equal

from pygrackle import \
    chemistry_data, \
    setup_fluid_container

from pygrackle.utilities.primordial_equilibrium import \
    total_cooling, \
    nHI, nHII, nHeI, nHeII, nHeIII

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs

def test_equilibrium():
    """
    Test ionization equilibrium for a primordial gas against the analytical
    solution using the same rates.

    This test takes a long time because the iteration to convergence near
    T = 1e4 K is slow.
    """

    my_chem = chemistry_data()
    my_chem.use_grackle = 1
    my_chem.with_radiative_cooling = 0
    my_chem.primordial_chemistry = 1
    my_chem.metal_cooling = 0
    my_chem.UVbackground = 0
    my_chem.comoving_coordinates = 0
    my_chem.a_units = 1.0
    my_chem.a_value = 1.0
    my_chem.density_units = mass_hydrogen_cgs
    my_chem.length_units = 1.0
    my_chem.time_units = 1.0
    my_chem.velocity_units = my_chem.length_units / my_chem.time_units
    fc = setup_fluid_container(my_chem, temperature=np.logspace(4.5, 9, 200),
                               converge=True, tolerance=1e-6, max_iterations=np.inf)

    fc.calculate_temperature()
    fc.calculate_cooling_time()
    t_sort = np.argsort(fc["temperature"])
    t_cool = fc["cooling_time"][t_sort] * my_chem.time_units
    my_T = fc["temperature"][t_sort]
    my_nH = fc.calculate_hydrogen_number_density().mean()

    cooling_rate_eq = -1*total_cooling(my_T, my_nH) / my_nH**2
    cooling_rate_g = fc["energy"][t_sort] / t_cool * fc["density"] * \
        my_chem.density_units / my_nH**2

    # Make a cooling rate figure
    fontsize = 14
    axes = pyplot.axes()
    axes.loglog(my_T, -1*cooling_rate_eq, color='black', alpha=0.7,
                linestyle="--", linewidth=1.5)
    axes.loglog(my_T, -1*cooling_rate_g, color='black', alpha=0.7,
                linestyle="-", linewidth=1)
    axes.xaxis.set_label_text('T [K]', fontsize=fontsize)
    axes.yaxis.set_label_text(r'$-\Lambda$ / n${_{\rm H}}^{2}$ [erg s$^{-1}$ cm$^{3}$]',
                              fontsize=fontsize)
    axes.set_xlim(1e4, 1e9)
    axes.set_ylim(1e-26, 2e-22)
    tick_labels = axes.xaxis.get_ticklabels() + \
        axes.yaxis.get_ticklabels()
    for tick_label in tick_labels:
        tick_label.set_size(fontsize)
    pyplot.savefig('cooling.png')
    pyplot.clf()

    # Make an ionization balance figure
    nH_eq = nHI(my_T, my_nH) + nHII(my_T, my_nH)
    nH_g  = fc["HI"] + fc["HII"]
    fHI_eq = nHI(my_T, my_nH) / nH_eq
    fHI_g  = fc["HI"] / nH_g
    fHII_eq = nHII(my_T, my_nH) / nH_eq
    fHII_g  = fc["HII"] / nH_g

    nHe_eq = nHeI(my_T, my_nH) + nHeII(my_T, my_nH) + nHeIII(my_T, my_nH)
    nHe_g  = fc["HeI"] + fc["HeII"] + fc["HeIII"]
    fHeI_eq = nHeI(my_T, my_nH) / nHe_eq
    fHeI_g  = fc["HeI"] / nHe_g
    fHeII_eq = nHeII(my_T, my_nH) / nHe_eq
    fHeII_g  = fc["HeII"] / nHe_g
    fHeIII_eq = nHeIII(my_T, my_nH) / nHe_eq
    fHeIII_g  = fc["HeIII"] / nHe_g

    # Plot H ions
    axes = pyplot.axes()
    axes.loglog(my_T, fHI_eq,  color="#B82E00", alpha=0.7, linestyle="--", linewidth=1.5)
    axes.loglog(my_T, fHII_eq, color="#B88A00", alpha=0.7, linestyle="--", linewidth=1.5)
    axes.loglog(my_T, fHI_g,  label="HI", color="#B82E00", alpha=0.7, linestyle="-", linewidth=1.)
    axes.loglog(my_T, fHII_g, label="HII", color="#B88A00", alpha=0.7, linestyle="-", linewidth=1.)

    # Plot He ions
    axes.loglog(my_T, fHeI_eq,   color="#002EB8", alpha=0.7, linestyle="--", linewidth=1.5)
    axes.loglog(my_T, fHeII_eq,  color="#008AB8", alpha=0.7, linestyle="--", linewidth=1.5)
    axes.loglog(my_T, fHeIII_eq, color="#00B88A", alpha=0.7, linestyle="--", linewidth=1.5)
    axes.loglog(my_T, fHeI_g,   label="HeI", color="#002EB8", alpha=0.7, linestyle="-", linewidth=1.)
    axes.loglog(my_T, fHeII_g,  label="HeII", color="#008AB8", alpha=0.7, linestyle="-", linewidth=1.)
    axes.loglog(my_T, fHeIII_g, label="HeIII", color="#00B88A", alpha=0.7, linestyle="-", linewidth=1.)

    axes.xaxis.set_label_text('T [K]', fontsize=fontsize)
    axes.yaxis.set_label_text('fraction', fontsize=fontsize)
    axes.set_xlim(1e4, 1e9)
    axes.set_ylim(1e-10, 1)
    tick_labels = axes.xaxis.get_ticklabels() + \
        axes.yaxis.get_ticklabels()
    for tick_label in tick_labels:
        tick_label.set_size(fontsize)
    axes.legend(loc='best', prop=dict(size=fontsize))
    pyplot.savefig('fractions.png')

    test_precision = 1

    # test the cooling rates
    cool_ratio = cooling_rate_eq / cooling_rate_g
    assert_rel_equal(cooling_rate_eq, cooling_rate_g, test_precision,
                     "Equilibrium cooling rates disagree with min/max = %f/%f." %
                     (cool_ratio.min(), cool_ratio.max()))

    # test the ionization balance
    assert_rel_equal(fHI_eq, fHI_g, test_precision,
                     "HI fractions disagree with min/max = %f/%f." %
                     ((fHI_eq / fHI_g).min(), (fHI_eq / fHI_g).max()))

    assert_rel_equal(fHII_eq, fHII_g, test_precision,
                     "HII fractions disagree with min/max = %f/%f." %
                     ((fHII_eq / fHII_g).min(), (fHII_eq / fHII_g).max()))

    assert_rel_equal(fHeI_eq, fHeI_g, test_precision,
                     "HeI fractions disagree with min/max = %f/%f." %
                     ((fHeI_eq / fHeI_g).min(), (fHeI_eq / fHeI_g).max()))

    assert_rel_equal(fHeII_eq, fHeII_g, test_precision,
                     "HeII fractions disagree with min/max = %f/%f." %
                     ((fHeII_eq / fHeII_g).min(), (fHeII_eq / fHeII_g).max()))

    assert_rel_equal(fHeIII_eq, fHeIII_g, test_precision,
                     "HeIII fractions disagree with min/max = %f/%f." %
                     ((fHeIII_eq / fHeIII_g).min(), (fHeIII_eq / fHeIII_g).max()))
