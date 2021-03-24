########################################################################
#
# Tests rates calculated when chemistry is initialized.
#
# Written by Ewan Jones, 19/03/2021
#
########################################################################

import h5py
import numpy as np 
import matplotlib.pyplot as plt 
import os

#* Import necessary functions from grackle.
from pygrackle import chemistry_data, setup_fluid_container

#* Import necessary constants from grackle.
from pygrackle.utilities.physical_constants import mass_hydrogen_cgs

#def test_rate_initialisation():


def test_rate_initialisation():
    """
    Test that the rate tables are initialized correctly.

    If writeHDF5 is set to true, all rate coefficients calculated in the initialisation routine will 
    be saved into a hdf5 file. This should be used for testing and is not meant for frequent use.
    """
    #* This should be False in all use cases. Set to True only for testing.
    writeHDF5 = False

    #* Initialise chemistry_data instance
    my_chemistry = chemistry_data()

    #* Set parameters
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 0
    my_chemistry.primordial_chemistry = 1
    my_chemistry.metal_cooling = 0
    my_chemistry.UVbackground = 0
    my_chemistry.comoving_coordinates = 0
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0
    my_chemistry.density_units = mass_hydrogen_cgs
    my_chemistry.length_units = 1.0
    my_chemistry.time_units = 1.0
    my_chemistry.velocity_units = my_chemistry.length_units / my_chemistry.time_units

    #* Create fluid container from the chemistry_data instance above. 
    #* This initialises the rate coefficients.
    my_fluidContainer = setup_fluid_container(my_chemistry, temperature=np.logspace(4.5, 9, 200),
                        converge=True, tolerance=1e-6, max_iterations=np.inf)
    
    #* Dictionary of all rate variables which will be checked.
    testRates = {"k1": my_chemistry.k1, "k2": my_chemistry.k2, "k3": my_chemistry.k3, "k4": my_chemistry.k4, "k5": my_chemistry.k5,
                "k6": my_chemistry.k6, "k7": my_chemistry.k7, "k8": my_chemistry.k8, "k9": my_chemistry.k9, "k10": my_chemistry.k10,
                "k11": my_chemistry.k1, "k12": my_chemistry.k12, "k13": my_chemistry.k13, "k14": my_chemistry.k14, "k15": my_chemistry.k15,
                "k16": my_chemistry.k16,"k17": my_chemistry.k17, "k18": my_chemistry.k18, "k19": my_chemistry.k19, "k20": my_chemistry.k20,
                "k21": my_chemistry.k21, "k22": my_chemistry.k22, "k23": my_chemistry.k23, "k24": my_chemistry.k24, "k25": my_chemistry.k25,
                "k26": my_chemistry.k26, "k27": my_chemistry.k27, "k28": my_chemistry.k28, "k29": my_chemistry.k29, "k30": my_chemistry.k30,
                "k31": my_chemistry.k31, "k50": my_chemistry.k50, "k51": my_chemistry.k51, "k52": my_chemistry.k52, "k53": my_chemistry.k53,
                "k54": my_chemistry.k54,"k55": my_chemistry.k55, "k56": my_chemistry.k56, "k57": my_chemistry.k57, "k58": my_chemistry.k58,
                "n_cr_n": my_chemistry.n_cr_n, "n_cr_d1": my_chemistry.n_cr_d1, "n_cr_d2": my_chemistry.n_cr_d2, "ceHI": my_chemistry.ceHI,
                "ceHeI": my_chemistry.ceHeI, "ceHeII": my_chemistry.ceHeII, "ciHI": my_chemistry.ciHI, "ciHeI": my_chemistry.ciHeI,
                "ciHeIS": my_chemistry.ciHeIS, "ciHeII": my_chemistry.ciHeII, "reHII": my_chemistry.reHII,"reHeII1": my_chemistry.reHeII1,
                "reHeII2": my_chemistry.reHeII2,"reHeIII": my_chemistry.reHeIII, "brem": my_chemistry.brem, "hyd01k": my_chemistry.hyd01k,
                "h2k01": my_chemistry.h2k01, "vibh": my_chemistry.vibh, "orth": my_chemistry.roth, "rotl": my_chemistry.rotl,
                "HDlte": my_chemistry.HDlte, "HDlow": my_chemistry.HDlow,"cieco": my_chemistry.cieco, "GAHI": my_chemistry.GAHI,
                "GAH2": my_chemistry.GAH2, "GAHe": my_chemistry.GAHe, "GAHp": my_chemistry.GAHp, "GAel": my_chemistry.GAel,
                "H2LTE": my_chemistry.H2LTE, "k13dd": my_chemistry.k13dd, "h2dust": my_chemistry.h2dust}

    #* Write rate coefficients to a hdf5 file if user specified such.
    if writeHDF5:
        #If file already exists delete it so that the new one can be created.
        if os.path.exists("initialised_rates.hdf5"):
            os.remove("initialised_rates.hdf5")

        #Write the file.
        f = h5py.File("initialised_rates.hdf5", "w-")
        for rate_key in testRates:
            f.create_dataset(rate_key, data=testRates[rate_key])
        f.close()

    #* Compare rates with the correct ones which are stored and check they are in agreement.
    passTest = True
    correctRates = h5py.File("tests/example_answers/correct_rates.hdf5", "r")
    for rate_key in testRates:
        if np.allclose(correctRates[rate_key], testRates[rate_key], atol=1e-10) == False:
            print("{} does not agree with correct rates.".format(rate_key))
            passTest = False
    
    assert(passTest == True)

    

        



    



    

            


