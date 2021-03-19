//Definition of function which saves all rates calculated in initialize_chemistry_data.c to various .txt files.
//This is used for testing purposes.
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

//Function definition.
int writeRates(char language[50], chemistry_data *my_chemistry, chemistry_data_storage *my_rates) {
    FILE *fp;

    //* Set up necessary folders and naming conventions.
    char directory[500] = "./ratesComparison/";
    char innerFolderName[20] = "_results/";
    strcat(language, innerFolderName);
    strcat(directory, language);

    //* Write k1-k58 rates.
    char fileName1[1000] = " ";
    strcpy(fileName1, directory);
    strcat(fileName1, "k1-k58_rates.txt");
    fp = fopen(fileName1, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->k1[line]); //0
        fprintf(fp, "%.12e,", my_rates->k2[line]); //1
        fprintf(fp, "%.12e,", my_rates->k3[line]); //2
        fprintf(fp, "%.12e,", my_rates->k4[line]); //3
        fprintf(fp, "%.12e,", my_rates->k5[line]); //4
        fprintf(fp, "%.12e,", my_rates->k6[line]); //5
        fprintf(fp, "%.12e,", my_rates->k7[line]); //6
        fprintf(fp, "%.12e,", my_rates->k8[line]); //7
        fprintf(fp, "%.12e,", my_rates->k9[line]); //8
        fprintf(fp, "%.12e,", my_rates->k10[line]); //9
        fprintf(fp, "%.12e,", my_rates->k11[line]); //10
        fprintf(fp, "%.12e,", my_rates->k12[line]); //11
        fprintf(fp, "%.12e,", my_rates->k13[line]); //12
        fprintf(fp, "%.12e,", my_rates->k14[line]); //13
        fprintf(fp, "%.12e,", my_rates->k15[line]); //14
        fprintf(fp, "%.12e,", my_rates->k16[line]); //15
        fprintf(fp, "%.12e,", my_rates->k17[line]); //16
        fprintf(fp, "%.12e,", my_rates->k18[line]); //17
        fprintf(fp, "%.12e,", my_rates->k19[line]); //18
        fprintf(fp, "%.12e,", my_rates->k20[line]); //19
        fprintf(fp, "%.12e,", my_rates->k21[line]); //20
        fprintf(fp, "%.12e,", my_rates->k22[line]); //21
        fprintf(fp, "%.12e,", my_rates->k23[line]); //22
        fprintf(fp, "%.12e,", my_rates->k50[line]); //23
        fprintf(fp, "%.12e,", my_rates->k51[line]); //24
        fprintf(fp, "%.12e,", my_rates->k52[line]); //25
        fprintf(fp, "%.12e,", my_rates->k53[line]); //26
        fprintf(fp, "%.12e,", my_rates->k54[line]); //27
        fprintf(fp, "%.12e,", my_rates->k55[line]); //28
        fprintf(fp, "%.12e,", my_rates->k56[line]); //29
        fprintf(fp, "%.12e,", my_rates->k57[line]); //30
        fprintf(fp, "%.12e,", my_rates->k58[line]); //31
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write H2 formation heating terms.
    char fileName2[1000] = " ";
    strcpy(fileName2, directory);
    strcat(fileName2, "H2formHeating_rates.txt");
    fp = fopen(fileName2, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->n_cr_n[line]);
        fprintf(fp, "%.12e,", my_rates->n_cr_d1[line]);
        fprintf(fp, "%.12e,", my_rates->n_cr_d2[line]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write cooling and heating rates.
    char fileName3[1000] = " ";
    strcpy(fileName3, directory);
    strcat(fileName3, "coolingAndHeating_rates.txt");
    fp = fopen(fileName3, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->ceHI[line]); //0
        fprintf(fp, "%.12e,", my_rates->ceHeI[line]); //1
        fprintf(fp, "%.12e,", my_rates->ceHeII[line]); //2
        fprintf(fp, "%.12e,", my_rates->ciHI[line]); //3
        fprintf(fp, "%.12e,", my_rates->ciHeI[line]); //4
        fprintf(fp, "%.12e,", my_rates->ciHI[line]); //5
        fprintf(fp, "%.12e,", my_rates->ciHeIS[line]); //6
        fprintf(fp, "%.12e,", my_rates->ciHeII[line]); //7
        fprintf(fp, "%.12e,", my_rates->reHII[line]); //8
        fprintf(fp, "%.12e,", my_rates->reHeII1[line]); //9
        fprintf(fp, "%.12e,", my_rates->reHeII2[line]); //10
        fprintf(fp, "%.12e,", my_rates->reHeIII[line]); //11
        fprintf(fp, "%.12e,", my_rates->brem[line]); //12
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write molecular hydrogen cooling rates.
    char fileName4[1000] = " ";
    strcpy(fileName4, directory);
    strcat(fileName4, "molecHydrogenCooling_rates.txt");
    fp = fopen(fileName4, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->hyd01k[line]);
        fprintf(fp, "%.12e,", my_rates->h2k01[line]);
        fprintf(fp, "%.12e,", my_rates->vibh[line]);
        fprintf(fp, "%.12e,", my_rates->roth[line]);
        fprintf(fp, "%.12e,", my_rates->rotl[line]);
        fprintf(fp, "%.12e,", my_rates->HDlte[line]);
        fprintf(fp, "%.12e,", my_rates->HDlow[line]);
        fprintf(fp, "%.12e,", my_rates->cieco[line]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write low density rates.
    char fileName5[1000] = " ";
    strcpy(fileName5, directory);
    strcat(fileName5, "lowDensity_rates.txt");
    fp = fopen(fileName5, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->GAHI[line]); //0
        fprintf(fp, "%.12e,", my_rates->GAH2[line]); //1
        fprintf(fp, "%.12e,", my_rates->GAHe[line]); //2
        fprintf(fp, "%.12e,", my_rates->GAHp[line]); //3
        fprintf(fp, "%.12e,", my_rates->GAel[line]); //4
        fprintf(fp, "%.12e,", my_rates->H2LTE[line]); //5
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write k13dd.
    char fileName6[1000] = " ";
    strcpy(fileName6, directory);
    strcat(fileName6, "k13dd.txt");
    fp = fopen(fileName6, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        for (int col = 0; col < 14; col++) {
            fprintf(fp, "%.12e,", my_rates->k13dd[line + col*my_chemistry->NumberOfTemperatureBins]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write h2dust.
    char fileName7[1000] = " ";
    strcpy(fileName7, directory);
    strcat(fileName7, "h2dust.txt");
    fp = fopen(fileName7, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        for (int col = 0; col < my_chemistry->NumberOfDustTemperatureBins; col++) {
            fprintf(fp, "%.12e,", my_rates->h2dust[line + col*my_chemistry->NumberOfTemperatureBins]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return SUCCESS;
}