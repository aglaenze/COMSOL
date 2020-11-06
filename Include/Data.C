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

//____________________________________________
void PrintList(
                 TString name,
                 Double_t* list,
                 Int_t num,
                 TString format = "%f" )
{
    //std::cout << "const Double_t " << name << "["<<size<< "] = {";
    std::cout << "Double_t " << name << "[" << num << "]" << " = {";
    for( Int_t i=0; i<num; ++i )
    {
        std::cout << Form( format.Data(), list[i] );
        if( i != num-1 ) std::cout << ", ";
        else std::cout << "};";
    }
    std::cout << std::endl;
}



//___________________________
// GAIN ZZBOT

Double_t gainListDataAriC4H10[11] = {430.621, 601.575, 835.974, 1167.004, 1626.864, 2250.949, 3148.231, 4456.768, 6222.110, 8731.924, 12189.061};
Double_t gainErrorListDataAriC4H10[11] = {13.817, 12.924, 11.392, 12.381, 13.891, 23.617, 26.841, 47.777, 57.497, 118.237, 140.382};
Double_t hvMeshListDataAriC4H10[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
Double_t hvDriftListDataAriC4H10[11] = {540.000, 550.000, 560.000, 570.000, 580.000, 590.000, 600.000, 610.000, 620.000, 630.000, 640.000};
Double_t resolutionListDataAriC4H10[11] = {35.433, 30.824, 28.267, 27.155, 26.954, 26.493, 26.466, 26.976, 27.315, 28.186, 29.380};
Double_t resolutionErrorListDataAriC4H10[11] = {12.108, 7.574, 5.653, 4.522, 3.770, 4.498, 3.839, 4.628, 3.993, 5.553, 4.886};

Double_t hvMeshListDataAriC4H10_old[11] = {340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440};
Double_t gainListDataAriC4H10_old[11] = {1098.851, 1510.038, 2128.435, 3014.917, 4277.629, 6067.106, 8638.780, 12246.708, 17504.688, 25397.665, 36880.590};
Double_t gainErrorListDataAriC4H10_old[11] = {0.308, 0.301, 0.501, 0.575, 0.782, 0.824, 1.175, 2.074, 2.602, 5.148, 7.196};



void LoadGainData(std::string gasName, const Int_t dataNum, Double_t hvMeshListData[], Double_t hvDriftListData[], Double_t gainListData[], Double_t gainErrorListData[], Double_t gainListData_old[], Double_t gainErrorListData_old[]) {

    if (gasName=="Ar-iC4H10") {
        for (int i = 0; i < dataNum; i++) {
            hvMeshListData[i] = hvMeshListDataAriC4H10[i];
            hvDriftListData[i] = hvDriftListDataAriC4H10[i];
            gainListData[i] = gainListDataAriC4H10[i];
            gainErrorListData[i] = gainErrorListDataAriC4H10[i];
            gainListData_old[i] = gainListDataAriC4H10_old[i];
            gainErrorListData_old[i] = gainErrorListDataAriC4H10_old[i];
        }
    }
    else {std::cout << "No data for this gas"; return;}
}



// ________________________
// IBF ZZBOT

Double_t hvMeshListIBF1AriC4H10[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
Double_t hvDriftListIBF1AriC4H10[11] = {540.000, 550.000, 560.000, 570.000, 580.000, 590.000, 600.000, 610.000, 620.000, 630.000, 640.000};
Double_t ibfRawAriC4H10[11] = {1.219, 1.818, 1.384, 0.978, 1.039, 1.005, 1.020, 0.949, 0.849, 0.723, 0.734};
Double_t ibfErrorAriC4H10[11] = {0.481, 0.328, 0.215, 0.125, 0.107, 0.096, 0.077, 0.041, 0.029, 0.025, 0.016};
Double_t ionBackFlowCorrectedVect1AriC4H10[11] = {1.216678, 1.816338, 1.382804, 0.977143, 1.038385, 1.004556, 1.019682, 0.948776, 0.848839, 0.722885, 0.733918};
Double_t ionBackFlowCorrectedErrorVect1AriC4H10[11] = {0.481075, 0.328036, 0.215016, 0.125009, 0.107005, 0.096005, 0.077003, 0.041002, 0.029001, 0.025002, 0.016001};

Double_t hvMeshListIBF1AriC4H10_old[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
Double_t ionBackFlowCorrectedVect1AriC4H10_old[11] = {0.881, 0.936, 0.888, 0.886, 0.876, 0.905, 0.876, 0.860, 0.832, 0.794, 0.749};
Double_t ionBackFlowCorrectedErrorVect1AriC4H10_old[11] = {0.076, 0.054, 0.038, 0.027, 0.019, 0.014, 0.009, 0.007, 0.005, 0.003, 0.003};


/*
Double_t ionBackFlowCorrectedVect1AriC4H10[11], ionBackFlowCorrectedErrorVect1AriC4H10[11];

void ComputeCorrectedIbfAr(Double_t ionBackFlowCorrectedVect1AriC4H10[], Double_t ionBackFlowCorrectedErrorVect1AriC4H10[]) {
    for (int i=0; i<11; i++) {
        ionBackFlowCorrectedVect1AriC4H10[i] = ibfRawAriC4H10[i] - 1./gainListDataAriC4H10[i];
        ionBackFlowCorrectedErrorVect1AriC4H10[i] = ibfErrorAriC4H10[i] + gainErrorListDataAriC4H10[i]/(gainListDataAriC4H10[i]*gainListDataAriC4H10[i]);
    }
    PrintList("ionBackFlowCorrectedVect1AriC4H10", ionBackFlowCorrectedVect1AriC4H10, 11);
    PrintList("ionBackFlowCorrectedErrorVect1AriC4H10", ionBackFlowCorrectedErrorVect1AriC4H10, 11);
    
}
*/


void LoadIbfData(std::string gasName, const Int_t dataNum, Double_t hvMeshListIBF1[], Double_t hvDriftListIBF1[], Double_t ionBackFlowCorrectedVect1[], Double_t ionBackFlowCorrectedErrorVect1[], Double_t ionBackFlowCorrectedVect1_old[], Double_t ionBackFlowCorrectedErrorVect1_old[]) {

    if (gasName=="Ar-iC4H10") {
        //ComputeCorrectedIbfAr(ionBackFlowCorrectedVect1, ionBackFlowCorrectedErrorVect1);
        for (int i = 0; i < dataNum; i++) {
            hvMeshListIBF1[i] = hvMeshListIBF1AriC4H10[i];
            hvDriftListIBF1[i] = hvDriftListIBF1AriC4H10[i];
            ionBackFlowCorrectedVect1[i] = ionBackFlowCorrectedVect1AriC4H10[i];
            ionBackFlowCorrectedErrorVect1[i] = ionBackFlowCorrectedErrorVect1AriC4H10[i];
            ionBackFlowCorrectedVect1_old[i] = ionBackFlowCorrectedVect1AriC4H10_old[i];
            ionBackFlowCorrectedErrorVect1_old[i] = ionBackFlowCorrectedErrorVect1AriC4H10_old[i];
        }
    }
    else {std::cout << "No data for this gas"; return;}
}

        
