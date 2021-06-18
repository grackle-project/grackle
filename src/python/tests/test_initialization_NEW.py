########################################################################
#
# Tests rates calculated when chemistry is initialized.
#
# Written by Ewan Jones, 19/03/2021
#
########################################################################

import h5py
import numpy as np
import os
from os.path import expanduser

#* Import necessary functions from grackle.
from pygrackle import chemistry_data, setup_fluid_container

#* Import necessary constants from grackle.
from pygrackle.utilities.physical_constants import mass_hydrogen_cgs
from python.tests.test_initialisation import test_rate_initialisation

#* Function which returns chemistry_data instance with default initialisation settings.
def get_defChem():
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

    return my_chemistry


#* Function which sets parameters for a given parameter set.
def set_parameters(parSet, my_chemistry):
    #Default parameter set.
    my_chemistry.CaseBRecombination = 0
    my_chemistry.h2_charge_exchange_rate = 1
    my_chemistry.h2dust_rate = 1
    my_chemistry.collisional_excitation_rates = 1
    my_chemistry.collisional_ionisation_rates = 1
    my_chemistry.recombination_cooling_rates = 1
    my_chemistry.bremsstrahlung_cooling_rates = 1
    my_chemistry.h2_h_cooling_rate == 1
    my_chemistry.photoelectric_heating = -1
    if parSet == 1:
        my_chemistry.three_body_rate = 0
        return True
    #Alternate parameter set.
    elif parSet == 2:
        my_chemistry.CaseBRecombination = 1
        my_chemistry.h2_charge_exchange_rate = 0
        my_chemistry.three_body_rate = 1
        my_chemistry.h2dust_rate = 0
        my_chemistry.collisional_excitation_rates = 0
        my_chemistry.collisional_ionisation_rates = 0
        my_chemistry.recombination_cooling_rates = 0
        my_chemistry.bremsstrahlung_cooling_rates = 0
        my_chemistry.h2_h_cooling_rate == 0
        my_chemistry.photoelectric_heating = 2
        return True
    #Default parameter sets with other three_body_rate values.
    elif parSet == 3:
        my_chemistry.three_body_rate = 2
        return True
    elif parSet == 4:
        my_chemistry.three_body_rate = 3
        return True
    elif parSet == 5:
        my_chemistry.three_body_rate = 4
        return True
    elif parSet == 6:
        my_chemistry.three_body_rate = 5
        return True
    #Invalid parameter set.
    else:
        return False


#* Function which prints the values of the parameter set in use.
def print_parameter_set(my_chemistry):
    parameters = ["CaseBRecombination", "h2_charge_exchange_rate", "three_body_rate", "h2dust_rate", "collisional_excitation_rates",\
                    "collisional_ionisation_rates", "recombination_cooling_rates", "bremsstrahlung_cooling_rates",\
                    "h2_h_cooling_rate", "photoelectric_heating"]
    for parameter in parameters:
        print(parameter + ":", getattr(my_chemistry, parameter))


#* Function which writes initialised rates for each parameter set to a hdf5 file.
def write_init_rates(rateNames):
    """
    Write initialized rates to a hdf5 file.
    """

    #* Navigate to the directory where the file is located.
    filePath = os.path.abspath(__file__)
    dirPath = os.path.dirname(filePath)
    os.chdir(dirPath)

    #* Calculate rates for each parameter set and write to hdf5 file
    #Create and open file. If the file already exists this will overwrite it.
    f = h5py.File("initialised_rates.h5", "w")
    for parSet in [1,2,3,4,5,6]:
        #Set chemistry parameters
        my_chem = get_defChem()
        if not set_parameters(parSet, my_chem):
            print("Parameter set could not be initialized.")
            exit()
        #Initialise and save rate coefficients.
        my_fluidContainer = setup_fluid_container(my_chem, temperature=np.logspace(4.5, 9, 200),
                            converge=True, tolerance=1e-6, max_iterations=np.inf)
        for rate in rateNames:
            f.create_dataset(rate + f"_{parSet}", data=getattr(my_chem, rate))
    #Close the file.
    f.close()

#* Function which checks that rates have been initialised correctly.
def test_initialization(rateNames, inputFile="initialised_rates.h5", correctFile="example_answers/correct_rates.h5", atol=0, rtol=1e-7):
    #* Open both files.
    testFile = h5py.File(inputFile, "r")
    correctFile = h5py.File(correctFile, "r")

    #* Check each rate
    incorrectRates = open("incorrectRates.txt", "w+")
    for parSet in [1,2,3,4,5,6]:
        for rateName in rateNames:
            rate = rateNames + f"_{parSet}"
            testRates = testFile[rate]
            correctRates = correctFile[rate]
            rateDiscrepancy = False
            for i, testVal in enumerate(testRates):
                if testVal != correctRates[i]:
                    print(rate + f"INCORRECT AT INDEX {i}")
                    rateDiscrepancy = True
                if rateDiscrepancy:
                    incorrectRates.write(rate + "\n")
    
    #Close the files
    incorrectRates.close()
    testFile.close()
    correctFile.close()
        

#* List of all rates that need to be checked.
rates = "k1,k3,k4,k2,k5,k6,k7,k8,k9,k10,k11,k12,k14,k15,k16,k17,k18,k19,k20,k23,"\
        "k13dd,k13,k21,k22,k50,k51,k52,k53,k54,k55,k56,k57,k58,h2dust,n_cr_n,n_cr_d1,"\
        "n_cr_d2,ceHI,ceHeI,ceHeII,ciHeIS,ciHI,ciHeI,ciHeII,reHII,reHeII1,reHeII2,"\
        "reHeIII,brem,vibh,hyd01k,h2k01,rotl,roth,GP99LowDensityLimit,GP99HighDensityLimit,"\
        "GAHI,GAH2,GAHe,GAHp,GAel,H2LTE,HDlte,HDlow,cieco,gas_grain,regr,comp,gammah,"\
        "gamma_isrf".split(',')

test_rate_initialisation(rates)