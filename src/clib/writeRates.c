//Definition of function which saves all rates calculated in initialize_chemistry_data.c to various .txt files.
//This is used for testing purposes.
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

//Function definition.
int writeRates(char language, chemistry_data *my_chemistry, chemistry_data_storage *my_rates) {
    FILE *fp;

    //* Set up necessary folders and naming conventions.
    char fileExtension[5] = ".txt";
    char directory[100] = "/ratesComparison/";
    char innerFolderName[50] = "_results";
    strcat(language, innerFolderName);
    strcat(directory, language);

    //* Write k1-k58 rates.
    char fileName1[100] = "k1-k58_rates_";
    strcat(fileName1, fileExtension);
    strcat(directory, fileName1);
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
    char fileName2[100] = "H2formHeating_rates_";
    strcat(fileName2, fileExtension);
    strcat(directory, fileName2);
    fp = fopen(fileName2, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->n_cr_n[line]);
        fprintf(fp, "%.12e,", my_rates->n_cr_d1[line]);
        fprintf(fp, "%.12e,", my_rates->n_cr_d2[line]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write cooling and heating rates.
    char fileName3[100] = "coolingAndHeating_rates_";
    strcat(fileName3, fileExtension);
    strcat(directory, fileName3);
    fp = fopen(fileName3, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->ceHI[line]);
        fprintf(fp, "%.12e,", my_rates->ceHeI[line]);
        fprintf(fp, "%.12e,", my_rates->ceHeII[line]);
        fprintf(fp, "%.12e,", my_rates->ciHI[line]);
        fprintf(fp, "%.12e,", my_rates->ciHeI[line]);
        fprintf(fp, "%.12e,", my_rates->ciHI[line]);
        fprintf(fp, "%.12e,", my_rates->ciHeIS[line]);
        fprintf(fp, "%.12e,", my_rates->ciHeII[line]);
        fprintf(fp, "%.12e,", my_rates->reHII[line]);
        fprintf(fp, "%.12e,", my_rates->reHeII1[line]);
        fprintf(fp, "%.12e,", my_rates->reHeII2[line]);
        fprintf(fp, "%.12e,", my_rates->reHeIII[line]);
        fprintf(fp, "%.12e,", my_rates->brem[line]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write molecular hydrogen cooling rates.
    char fileName4[100] = "molecHydrogenCooling_rates_";
    strcat(fileName4, fileExtension);
    strcat(directory, fileName4);
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
    char fileName5[100] = "lowDensity_rates_";
    strcat(fileName5, fileExtension);
    strcat(directory, fileName5);
    fp = fopen(fileName5, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        fprintf(fp, "%.12e,", my_rates->GAHI[line]);
        fprintf(fp, "%.12e,", my_rates->GAH2[line]);
        fprintf(fp, "%.12e,", my_rates->GAHe[line]);
        fprintf(fp, "%.12e,", my_rates->GAHp[line]);
        fprintf(fp, "%.12e,", my_rates->GAel[line]);
        fprintf(fp, "%.12e,", my_rates->H2LTE[line]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write k13dd.
    char fileName6[100] = "k13dd_";
    strcat(fileName6, fileExtension);
    strcat(directory, fileName6);
    fp = fopen(fileName6, "w");
    for (int line = 0; line < my_chemistry->NumberOfTemperatureBins; line++) {
        for (int col = 0; col < 14; col++) {
            fprintf(fp, "%.12e,", my_rates->k13dd[line + col*my_chemistry->NumberOfTemperatureBins]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    //*Write h2dust.
    char fileName7[100] = "h2dust_";
    strcat(fileName7, fileExtension);
    strcat(directory, fileName7);
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