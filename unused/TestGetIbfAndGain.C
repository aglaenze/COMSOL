#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "_Utils.C"



void TestGetIbfAndGain() {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    //____________________

    time_t t0 = time(NULL);
    gStyle->SetOptStat(0);
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
     
    // Get number of files to look at
    Int_t num = 1;
    //Int_t num = GetNumberOfFiles(path, "signal");


    TFile* fSignal = TFile::Open("rootFiles/Ar-iC4H10/model1/signal-340-540-2.root", "READ");
    
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    Int_t nEvents = tAvalanche->GetEntries();
    
    std::map <std::string, int> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    
    TTree* tSignalPad = (TTree*)fSignal->Get(Form("tSignal_%d", electrodeMap["pad"]));
    TTree* tSignalDrift = (TTree*)fSignal->Get(Form("tSignal_%d", electrodeMap["drift"]));
    TTree* tSignalMesh = (TTree*)fSignal->Get(Form("tSignal_%d", electrodeMap["mesh"]));

        
    int nEntries = tSignalPad->GetEntries();
    if (tSignalPad->GetEntries() != tSignalDrift->GetEntries() ) {
        std::cout << "not the same number of entries in tree pad and tree drift" << std::endl;
        return;
    }
    
    Double_t currentDriftIntMax = tSignalDrift->GetMaximum("totalCurrentInt");
    Double_t currentPadIntMin = tSignalPad->GetMinimum("totalCurrentInt");
    Double_t currentMeshIntMax = tSignalMesh->GetMaximum("totalCurrentInt");
    std::cout << "currentDriftIntMax = " << currentDriftIntMax << std::endl;
    std::cout << "currentMeshIntMax = " << currentMeshIntMax << std::endl;
    std::cout << "currentPadIntMin = " << currentPadIntMin << std::endl;
    std::cout << "IBF (reading the pads) = " << -currentDriftIntMax/currentPadIntMin*100 << "%" << std::endl;
    std::cout << "IBF (reading the mesh) = " << currentDriftIntMax/currentMeshIntMax*100 << "%" << std::endl;
    
    std::cout << std::endl;
    
    std::cout << "Number of charges in the drift = " << int(currentDriftIntMax/nEvents) << std::endl;
    std::cout << "Number of charges in the mesh = " << int(currentMeshIntMax/nEvents) << std::endl;
    std::cout << "Number of charges in the pads = " << -int(currentPadIntMin/nEvents) << std::endl;
    
    std::cout << std::endl;
    
    std::cout << "Number of electrons in the drift = " << int(tSignalDrift->GetMaximum("electronCurrentInt")/nEvents) << std::endl;
    std::cout << "Number of electrons in the mesh = " << int(tSignalMesh->GetMaximum("electronCurrentInt")/nEvents) << std::endl;
    std::cout << "Number of electrons in the pads = " << -int(tSignalPad->GetMinimum("electronCurrentInt")/nEvents) << std::endl;
    
    std::cout << std::endl;
    
    std::cout << "Number of ions in the drift = " << int(tSignalDrift->GetMaximum("ionCurrentInt")/nEvents) << std::endl;
    std::cout << "Number of ions in the mesh = " << int(tSignalMesh->GetMaximum("ionCurrentInt")/nEvents) << std::endl;
    std::cout << "Number of ions in the pads = " << -int(tSignalPad->GetMinimum("ionCurrentInt")/nEvents) << std::endl;

    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    return;
}


