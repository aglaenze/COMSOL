#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>

using namespace std;

int FindElectrodeIndexStart(int modelNum, double zPos) {
    
    vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    for (int i = 0; i < (int)zElectrodes.size()-1; i++) {
        if (zPos < zElectrodes[i] && zPos > zElectrodes[i+1]) return i;
    }
    return -1;
}

void DrawElectronsStart(const int modelNum = 15, TString fSignalName="", double zmin = 0, double zmax = 0) {
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    
    std::vector<float> *electronStartPointsInput = 0;
    tAvalanche->SetBranchAddress("electronStartPoints", &electronStartPointsInput);
    
    std::vector<float> electronStartPoints = {};
    TH1F* zElDistribution = new TH1F("hZelectrons", "Start z of electrons", 1000, zmin, zmax);
    const int nAvalanche = tAvalanche->GetEntries();
    for (int k = 0; k < nAvalanche; k++) {
        tAvalanche->GetEntry(k);
        electronStartPoints = *electronStartPointsInput;
        for (int j = 0; j< (int)electronStartPoints.size(); j++) {
            zElDistribution->Fill(electronStartPoints[j]);
        }
        electronStartPoints.clear();
    }
    //zElDistribution->Scale(1/zElDistribution->GetMaximum());
    //if (!(modelNum==1 || (modelNum >15 && modelNum < 19)) ) zElDistribution->SetMaximum(0.05* zElDistribution->GetMaximum());
    zElDistribution->GetXaxis()->SetTitle("z (cm)");
    //if (zmax < 0.1) zElDistribution->GetXaxis()->SetMaxDigits(2);
    
    zElDistribution->Draw("hist");
}

void DrawElectrodeStart(const int modelNum = 15, TString fSignalName="") {
    
    gStyle->SetLabelSize(.065, "X");
    gStyle->SetLabelOffset(0.01, "X");
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    
    int electrodeNum = GetElectrodeNum(modelNum);
    electrodeNum--; // on se fiche des pads
    
    map <string, int, NoSorting> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    string electrodeNames[electrodeNum];
    char *elNameChar[electrodeNum];
    map<string, int>::iterator it = electrodeMap.begin();
    int i = 0;
    while (it != electrodeMap.end()) {
        if (it->first != "pad") {
            electrodeNames[i] = it->first;
            elNameChar[i] = Form("%s", (it->first).c_str());
            i++;
        }
        it++;
    }
    
    int nStartElectrons[electrodeNum];
    for (int k = 0; k < electrodeNum; k++) {
        nStartElectrons[k] = 0;
    }
    
    std::vector<float> *electronStartPointsInput = 0;
    tAvalanche->SetBranchAddress("electronStartPoints", &electronStartPointsInput);
    
    std::vector<float> electronStartPoints = {};
    TH1F* elDistribution = new TH1F("hel", "Electrode below which electrons start", electrodeNum, 0, electrodeNum);
    const int nAvalanche = tAvalanche->GetEntries();
    for (int k = 0; k < nAvalanche; k++) {
        tAvalanche->GetEntry(k);
        electronStartPoints = *electronStartPointsInput;
        for (int j = 0; j< (int)electronStartPoints.size(); j++) {
            int elIdx = FindElectrodeIndexStart(modelNum, electronStartPoints[j]);
            if (elIdx < 0) cout << "problem with index < 0" << endl;
            nStartElectrons[elIdx]++;
        }
        electronStartPoints.clear();
    }
    for (int k = electrodeNum-1; k >= 0; k--) {
        elDistribution->Fill(elNameChar[k], (double)nStartElectrons[k]/nAvalanche);
    }
    
    // Combine GEM up and GEM down if there is a GEM
    string sGem = "GEM";
    int newDim = electrodeNum;
    bool gemPresent = false;
    for (int k = 0; k < electrodeNum; k++) {
        if (electrodeNames[k].find(sGem) != string::npos) {
            gemPresent = true;
            newDim = electrodeNum-1;
            break;
        }
    }
    
    // Change the dimension of elNameChar and nIonsStop to sum both GEM contributions
    char* elNameCharCopy[newDim];
    string electrodeNamesCopy[newDim];
    int nStartElectronsCopy[newDim];
    // initialisation of nIonsStopCopy et des histos
    for (int k = 0; k < newDim; k++) {
        nStartElectronsCopy[k] = 0;
    }
    
    // Loop to rename electrodes
    int idx1 = 0;
    for (int k = 0; k < electrodeNum; k++) {
        if (electrodeNames[k].find(sGem) != string::npos) {
            elNameCharCopy[idx1] = Form("%s", sGem.c_str());
            electrodeNamesCopy[idx1] = sGem;
        }
        else {
            if (k > 0 && electrodeNames[k-1].find(sGem) != string::npos) {idx1++;}
            elNameCharCopy[idx1] = elNameChar[k];
            electrodeNamesCopy[idx1] = electrodeNames[k];
            idx1++;
        }
    }
    
    // Loop to re-count electrons
    idx1 = 0;
    for (int k = 0; k < electrodeNum; k++) {
        if (k > 0 && electrodeNames[k].find(sGem)!= std::string::npos && electrodeNames[k-1].find(sGem) != std::string::npos) {idx1--;}
        nStartElectronsCopy[idx1] += nStartElectrons[k];
        idx1++;
    }
    
    TH1F* elDistributionCopy = new TH1F("hel2", "Electrode below which electrons start", newDim, 0, newDim);
    for (int k = newDim-1; k >= 0; k--) {
        elDistributionCopy->Fill(elNameCharCopy[k], (double)nStartElectronsCopy[k]/nAvalanche);
    }
    
    //elDistribution->Draw("hist");
    elDistributionCopy->Draw("hist");
}
