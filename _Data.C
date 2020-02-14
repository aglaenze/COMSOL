#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>


Int_t dataQuantity(std::string gasName) {
    if (gasName=="Ar-iC4H10") return 11;
    else {std::cout << "No data for this gas"; return 0;}
}


//___________________________
// GAIN ZZBOT

Double_t gainListDataAriC4H10[11] = {430.621, 601.575, 835.974, 1167.004, 1626.864, 2250.949, 3148.231, 4456.768, 6222.110, 8731.924, 12189.061};
Double_t gainErrorListDataAriC4H10[11] = {13.817, 12.924, 11.392, 12.381, 13.891, 23.617, 26.841, 47.777, 57.497, 118.237, 140.382};
Double_t hvMeshListDataAriC4H10[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
Double_t resolutionListDataAriC4H10[11] = {35.433, 30.824, 28.267, 27.155, 26.954, 26.493, 26.466, 26.976, 27.315, 28.186, 29.380};
Double_t resolutionErrorListDataAriC4H10[11] = {12.108, 7.574, 5.653, 4.522, 3.770, 4.498, 3.839, 4.628, 3.993, 5.553, 4.886};



void LoadGainData(std::string gasName, const Int_t dataNum, Double_t hvMeshListData[], Double_t gainListData[], Double_t gainErrorListData[]) {

    if (gasName=="Ar-iC4H10") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListData[i] = hvMeshListDataAriC4H10[i];
            gainListData[i] = gainListDataAriC4H10[i];
            gainErrorListData[i] = gainErrorListDataAriC4H10[i];
        }
    }
    else {std::cout << "No data for this gas"; return;}
}



// ________________________
// IBF ZZBOT

Double_t hvMeshListIBF1AriC4H10[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
Double_t ibfRawAriC4H10[11] = {1.219, 1.818, 1.384, 0.978, 1.039, 1.005, 1.020, 0.949, 0.849, 0.723, 0.734};
Double_t ibfErrorAriC4H10[11] = {0.481, 0.328, 0.215, 0.125, 0.107, 0.096, 0.077, 0.041, 0.029, 0.025, 0.016};

Double_t ionBackFlowCorrectedVect1AriC4H10[11], ionBackFlowCorrectedErrorVect1AriC4H10[11];

for (int i=0; i<11; i++) {
    ionBackFlowCorrectedVect1AriC4H10[i] = ibfRawAriC4H10[i] - 1./gainListDataAriC4H10[i];
    ionBackFlowCorrectedErrorVect1AriC4H10[i] = ibfErrorAriC4H10[i] + gainErrorListDataAriC4H10[i]/(gainListDataAriC4H10[i]*gainListDataAriC4H10[i]);
}



void LoadIbfData(std::string gasName, const Int_t dataNum, Double_t hvMeshListIBF1[], Double_t ionBackFlowCorrectedVect1[], Double_t ionBackFlowCorrectedErrorVect1[]) {

    if (gasName=="Ar-iC4H10") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListIBF1[i] = hvMeshListIBF1AriC4H10[i];
            ionBackFlowCorrectedVect1[i] = ionBackFlowCorrectedVect1AriC4H10[i];
            ionBackFlowCorrectedErrorVect1[i] = ionBackFlowCorrectedErrorVect1AriC4H10[i];
        }
    }
    else {std::cout << "No data for this gas"; return;}
}

        
