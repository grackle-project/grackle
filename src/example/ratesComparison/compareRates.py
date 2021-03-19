import datetime
import math 

#? Flags that can be set by the user to control terminal outputs.
showOOM = False
showUniqueRates = True
#?-----------------------------------------------

#Names of files which rates are stored in.
fileNames = ["k1-k58_rates.txt", "H2formHeating_rates.txt", "coolingAndHeating_rates.txt", "molecHydrogenCooling_rates.txt", "lowDensity_rates.txt", "k13dd.txt", "h2dust.txt"]

#File which detailed discrepancy data will be written to.
f = open("discrepanyLog.txt", "w+")

#Loop through all files.
for fileName in fileNames:
    #Write file header.
    f.write("Ouptuts compared on --> {} \n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    f.write("______________________________ \n")
    f.write("DISCREPANCIES FOUND IN {} \n".format(fileName.upper()))
    f.write("______________________________ \n")

    #Set file-specific variables to their initial values.
    discrepancy = 0 #Says if there is a discrepancy in the file.
    discrepancyMags = [] #List to hold all unique order-of-magnitudes in the discrepancies.
    discrepancyUniqueElements = {} #Dictionary which holds the element index of each discrepant rate, and attached to it all OOM at which it occurs

    #Open the C and Fortran rates.
    cFile_lines = open("c_results/{}".format(fileName), 'r').readlines()
    fortranFile_lines = open("fortran_results/{}".format(fileName), 'r').readlines()

    #Check that the files contain the same number of elements. If they do then go ahead with comparison.
    if (len(cFile_lines) != len(fortranFile_lines)):
        f.write("Differing number of lines \n")
    else:
        #Compare files line-by-line.
        for line in range(len(cFile_lines)):
            f.wroteLine = False #Flag which says if line declaration has been written in the file.
            #See if line as a whole differs between files. If it does then check element-by-element.
            if (cFile_lines[line] != fortranFile_lines[line]):
                discrepantRates = [] #Where elements will be stored if they are discrepant.

                #Get individual elements from line and compare them.
                cRates = cFile_lines[line].split(",")
                fortranRates = fortranFile_lines[line].split(",")
                for rateInd in range(len(cRates)):
                    #Account for edge cases.
                    if cRates[rateInd] == "\n" or (float(cRates[rateInd]) == 0. and float(fortranRates[rateInd]) == 0.):
                        break
                    #Check for discrepancy between two values.
                    else:
                        #Calculate relative discrepancy.
                        discrepancyVal = abs(float(cRates[rateInd]) - float(fortranRates[rateInd]))/max(float(fortranRates[rateInd]), float(cRates[rateInd]))

                        #If discrepancy is not deemed allowable do the following.
                        if discrepancyVal > 1e-10:
                            discrepancy = 1 #There is a discrepancy in the file.

                            #Calculate order-of-magnitude of the relative discrepancy and store them if they are unique.
                            OOM = "1e{}".format(math.floor(math.log10(discrepancyVal)))
                            if OOM not in discrepancyMags:
                                discrepancyMags.append(OOM)

                            #Store unique rate elements and all the OOM at which the discrepancy exists.
                            if rateInd not in discrepancyUniqueElements:
                                discrepancyUniqueElements[rateInd] = []
                            else:
                                if OOM not in discrepancyUniqueElements[rateInd]:
                                    discrepancyUniqueElements[rateInd].append(OOM)

                            #If line decleration has not been written then write it.
                            if f.wroteLine == False:
                                f.write("Line {}: \n".format(line))
                                f.writeedLine = True

                            #Add the element index of the discrepant value to the list.
                            discrepantRates.append(rateInd)

                            #Write the full details of the discrepancy to the output file.
                            f.write(" -Element {}. C value: {}. Fortran value {}. Discrepancy {} \n".format(rateInd, cRates[rateInd], fortranRates[rateInd], discrepancyVal))
    
    #Terminal outputs controlled by flags at top of file.

    #Shows orders of magnitude.
    if showOOM:
        print("Unique magnitudes of discrepancies in {} are: \n".format(fileName))
        print(discrepancyMags)
        print("----------------------------------------------------- \n")

    #Shows unique rates and their order of magnitude.
    if showUniqueRates:
        print("Elements that are discrepant in {} are: \n".format(fileName))
        for key in discrepancyUniqueElements:
            print("Element: {}      Magnitudes: {} \n".format(key, discrepancyUniqueElements[key], reverse=True))
        print("----------------------------------------------------- \n")
    

f.close()
                            
    
                
            

