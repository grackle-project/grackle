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

#* Import necessary functions from grackle.
from pygrackle import chemistry_data, setup_fluid_container

#* Import necessary constants from grackle.
from pygrackle.utilities.physical_constants import mass_hydrogen_cgs

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
    parameters = ["CaseBRecombination", "k11_rate", "three_body_rate", "h2dust_rate", "collisional_excitation_rates",\
                    "collisional_ionisation_rates", "recombination_cooling_rates", "bremsstrahlung_cooling_rates",\
                    "h2_h_cooling_rate", "photoelectric_heating"]
    for parameter in parameters:
        print(parameter + ":", getattr(my_chemistry, parameter))

#* Function which tells you unique order of magnitude discrepancies for a given rate between two hdf5 files.
def oom_discrepancies(rateName, parameterSet, rateFile1, rateFile2):
    name = rateName + f'_{parameterSet}'
    rates1 = rateFile1[name]
    rates2 = rateFile2[name]

    uniqueDiscrepancies = {}
    #Check each rate coefficient and calculate OOM discrepancy.
    for i in range(len(rates1)):
        if rates1[i] == 0 or rates2[i] == 0:
            None
        else:
            relDisc = abs(rates1[i] - rates2[i]) / rates1[i]
            if relDisc == 0:
                None
            else:
                DiscOOM = np.floor(np.log10(relDisc))
                if f'e^{DiscOOM}' not in uniqueDiscrepancies:
                    uniqueDiscrepancies[f'e^{DiscOOM}'] = []
                uniqueDiscrepancies[f'e^{DiscOOM}'].append(i)

    return uniqueDiscrepancies


#* Function which tests that the rates have been initialised correctly for each parameter set.
def test_rate_initialisation(printParameters=False, printOOMdiscrepanices=False):
    """
    Test that the rate tables are initialized correctly.

    If writeHDF5 is set to true, all rate coefficients calculated in the initialisation routine will 
    be saved into a hdf5 file. This should be used for testing and is not meant for frequent use.

    printParameters (bool) --> If set to True will print the parameter settings for each parameter set.
    """

    #* Navigate to the directory where the file is located.
    filePath = os.path.abspath(__file__)
    dirPath = os.path.dirname(filePath)
    os.chdir(dirPath)

    #* List of all rate variable names which will be checked.
    testRates = ["k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13",
                 "k14", "k15", "k16","k17", "k18", "k19", "k20", "k21", "k22", "k23", "k24", "k25",
                 "k26", "k27", "k28", "k29", "k30", "k31", "k50", "k51", "k52", "k53", "k54", "k55",
                 "k56", "k57", "k58", "n_cr_n", "n_cr_d1", "n_cr_d2", "ceHI", "ceHeI", "ceHeII",
                 "ciHI", "ciHeI", "ciHeIS", "ciHeII", "reHII", "reHeII1", "reHeII2","reHeIII", "brem",
                 "hyd01k", "h2k01", "vibh", "roth", "rotl", "HDlte", "HDlow", "cieco", "GAHI", "GAH2",
                 "GAHe", "GAHp", "GAel", "H2LTE", "k13dd", "h2dust"]

    #* Calculate rates for each parameter set and write to hdf5 file
    #Create and open file. If the file already exists this will overwrite it.
    f = h5py.File("initialised_rates.h5", "w")

    #Iterate over parameter sets.
    parSets = [1,2,3,4,5,6]
    for parSet in parSets:
        my_chemistry = get_defChem()
        #Set chemistry parameters.
        if set_parameters(parSet, my_chemistry):
            if printParameters:
                print(f"---Parameter Set {parSet}---\n")
                print_parameter_set(my_chemistry)
            #Initialise the rate coefficients.
            my_fluidContainer = setup_fluid_container(my_chemistry, temperature=np.logspace(4.5, 9, 200),
                                converge=True, tolerance=1e-6, max_iterations=np.inf)
        else:
            print("Invalid parameter set encountered.")
            exit()

        #Write rates to file.
        for rate_key in testRates:
            f.create_dataset(rate_key + f"_{parSet}", data=getattr(my_chemistry, rate_key))

    #Close the file.
    f.close()

    #* Compare rates with the correct ones which are stored and check they are in agreement
    correctRates = h5py.File("example_answers/correct_rates.h5", "r")
    initialisedRates = h5py.File("initialised_rates.h5", "r")
    
    #Print order-of-magnitude of discrepancies found.
    if printOOMdiscrepanices:
        print(oom_discrepancies('k13', 3, correctRates, initialisedRates))

    #Check all rates for each parameter set.
    for rate_key in testRates:
        for parSet in parSets:
            rate_name = rate_key + f"_{parSet}"
            #Check rates agree to what we deem is an acceptable relative tolerance.
            assert np.allclose(correctRates[rate_name], initialisedRates[rate_name], rtol=1e-7),\
                                f"Rate Coefficient {rate_name} does not agree. \n Correct rate:\
                                    {correctRates[rate_name][300]} \n Initialised rate: {initialisedRates[rate_name][300]} \n"