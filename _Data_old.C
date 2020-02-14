#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>


Int_t dataQuantity(std::string gasName) {
    if (gasName=="Ne") return 11;
    else if (gasName=="Ar-iC4H10") return 11;
    else if (gasName=="Ar-CO2") return 5;
    // else if (gasName=="Ar-CO2") return 4;    // Maxence's data
    else {std::cout << "Unknown gas"; return 0;}
}


//___________________________
// GAIN ZZBOT

// Ne
Double_t hvMeshListDataNe[11] = {340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440};
Double_t gainListDataNe[11] = {596.076, 751.592, 945.953, 1196.154, 1531.500, 1933.089, 2437.905, 3142.112, 4002.650, 5086.517, 6345.390};
Double_t gainErrorListDataNe[11] = {0.281, 0.180, 0.297, 0.302, 0.358, 0.398, 0.495, 0.606, 0.748, 0.916, 1.118};

//Ar-iC4H10
Double_t hvMeshListDataAriC4H10[11] = {340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440};
Double_t gainListDataAriC4H10[11] = {1098.851, 1510.038, 2128.435, 3014.917, 4277.629, 6067.106, 8638.780, 12246.708, 17504.688, 25397.665, 36880.590};
Double_t gainErrorListDataAriC4H10[11] = {0.308, 0.301, 0.501, 0.575, 0.782, 0.824, 1.175, 2.074, 2.602, 5.148, 7.196};

// GAIN LITTLE CHINESE
// Ar-CO2
/*
// Maxence's data
Double_t hvMeshListData[4] = {350, 360, 370, 380};
Double_t gainListData[4] = {1074.38, 1346.71, 1656.42, 2174.37};
Double_t gainErrorListData[4] = {0., 0., 0., 0.};
*/
Double_t gainListDataArCO2[5] = {1204.974, 2444.215, 5041.168, 10563.040, 22966.871};
Double_t gainErrorListDataArCO2[5] = {13.391, 18.834, 40.020, 51.919, 111.320};
//hvMeshDownList[5] = {350.000, 360.000, 370.000, 380.000, 390.000};
Double_t hvMeshListDataArCO2[5] = {350.000, 360.000, 370.000, 380.000, 390.000};
//hvMeshUpList[5] = {1050.000, 1080.000, 1110.000, 1140.000, 1170.000};


void LoadGainData(std::string gasName, const Int_t dataNum, Double_t hvMeshListData[], Double_t gainListData[], Double_t gainErrorListData[]) {

    if (gasName=="Ne") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListData[i] = hvMeshListDataNe[i];
            gainListData[i] = gainListDataNe[i];
            gainErrorListData[i] = gainErrorListDataNe[i];
        }
    }
    else if (gasName=="Ar-iC4H10") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListData[i] = hvMeshListDataAriC4H10[i];
            gainListData[i] = gainListDataAriC4H10[i];
            gainErrorListData[i] = gainErrorListDataAriC4H10[i];
        }
    }
    else if (gasName=="Ar-CO2") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListData[i] = hvMeshListDataArCO2[i];
            gainListData[i] = gainListDataArCO2[i];
            gainErrorListData[i] = gainErrorListDataArCO2[i];
        }
    }
    else {std::cout << "Unknown gas"; return;}
}



// ________________________
// IBF


// ZZBOT
// Ne
 Double_t hvMeshListIBF1Ne[11]  = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
 Double_t ionBackFlowCorrectedVect1Ne[11]  = {0.885, 0.858, 0.886, 0.900, 0.914, 0.898, 0.887, 0.866, 0.843, 0.812, 0.773};
 Double_t ionBackFlowCorrectedErrorVect1Ne[11]  = {0.048, 0.034, 0.025, 0.017, 0.012, 0.009, 0.006, 0.004, 0.003, 0.002, 0.002};

// Ar-iC4H10
Double_t hvMeshListIBF1AriC4H10[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
Double_t ionBackFlowCorrectedVect1AriC4H10[11] = {0.881, 0.936, 0.888, 0.886, 0.876, 0.905, 0.876, 0.860, 0.832, 0.794, 0.749};
Double_t ionBackFlowCorrectedErrorVect1AriC4H10[11] = {0.076, 0.054, 0.038, 0.027, 0.019, 0.014, 0.009, 0.007, 0.005, 0.003, 0.003};

// LITTLE CHINESE
// Ar-CO2
Double_t hvMeshListIBF1ArCO2[2] = {0., 0.};
Double_t ionBackFlowCorrectedVect1ArCO2[2] = {0., 0.};
Double_t ionBackFlowCorrectedErrorVect1ArCO2[2] = {0., 0.};

void LoadIbfData(std::string gasName, const Int_t dataNum, Double_t hvMeshListIBF1[], Double_t ionBackFlowCorrectedVect1[], Double_t ionBackFlowCorrectedErrorVect1[]) {

    if (gasName=="Ne") {
        for (int i = 0; i < dataNum; i++) {
                hvMeshListIBF1[i] = hvMeshListIBF1Ne[i];
                ionBackFlowCorrectedVect1[i] = ionBackFlowCorrectedVect1Ne[i];
                ionBackFlowCorrectedErrorVect1[i] = ionBackFlowCorrectedErrorVect1Ne[i];
            }
        }
    else if (gasName=="Ar-iC4H10") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListIBF1[i] = hvMeshListIBF1AriC4H10[i];
            ionBackFlowCorrectedVect1[i] = ionBackFlowCorrectedVect1AriC4H10[i];
            ionBackFlowCorrectedErrorVect1[i] = ionBackFlowCorrectedErrorVect1AriC4H10[i];
        }
    }
    else if (gasName=="Ar-CO2") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListIBF1[i] = hvMeshListIBF1ArCO2[i];
            ionBackFlowCorrectedVect1[i] = ionBackFlowCorrectedVect1ArCO2[i];
            ionBackFlowCorrectedErrorVect1[i] = ionBackFlowCorrectedErrorVect1ArCO2[i];
        }
    }
    else {std::cout << "Unknown gas"; return;}
}

        
