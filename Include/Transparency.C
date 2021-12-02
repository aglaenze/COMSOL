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
#include <TLine.h>
#include <TMath.h>

using namespace std;

void DrawLines(TH1F* hist) {
    
    int nBins = hist->GetNbinsX();
    double yMin = hist->GetMinimum();
    
    // first: draw the bottom horizontal line
    TLine* lh = new TLine(0, yMin, nBins, yMin);
    lh->Draw();
    
    // draw vertical lines
    TLine* lv[nBins+1];
    for (int i = 0; i < nBins+1; i++) {
        double yPos = yMin;
        if (i<nBins) yPos = hist->GetBinContent(i+1);
        else yPos = hist->GetBinContent(i);
        lv[i] = new TLine(i, yMin, i, yPos);
    }
    
    for (int i = 0; i < nBins+1; i++) {lv[i]->Draw();}
}

void DrawDyingIons(const int modelNum = 15, TString fSignalName=""){
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    
    /*
     gStyle->SetTitleFontSize(.06);
     gStyle->SetTitleSize(.06);
     
     gStyle->SetOptStat(0);
     gStyle->SetTitleFontSize(.05);
     gStyle->SetTitleXSize(.05);
     gStyle->SetTitleYSize(.05);
     gStyle->SetLabelSize(.05, "XY");
     gStyle->SetMarkerSize(0.3);
     gStyle->SetTextSize(0.05);
     */
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    map <string, int, NoSorting> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    int nAval = tAvalanche->GetEntries();
    
    int nElectrodes = electrodeMap.size()-2;    // on s'en fiche de la drift et des pads
    cout << "There are " << nElectrodes << " electrodes to look at" << endl;
    TH1F* hDyingIons[nElectrodes];
    
    // Initialisation des histos de tranparence
    int i = 0;
    string electrodeNames[nElectrodes];
    map<string, int>::iterator it = electrodeMap.begin();
    it++;    // ignore drift electrode
    while (it != electrodeMap.end()) {
        if (it->first != "pad") {
            electrodeNames[i] = it->first;
            hDyingIons[i] = new TH1F(Form("hDyingIonsFromBelow%s", electrodeNames[i].c_str()), "Dying ions", int(ddrift*1.1*200), 0, ddrift*1.1);
            i++;
        }
        it++;
    }
    
    vector<float> *ionStartPointsInput = 0, *ionEndPointsInput = 0;
    tAvalanche->SetBranchAddress("ionStartPoints", &ionStartPointsInput);
    tAvalanche->SetBranchAddress("ionEndPoints", &ionEndPointsInput);
    
    vector<float> ionStartPoints = {}, ionEndPoints = {};
    
    int nIons[nElectrodes];
    for (int k = 0; k < nElectrodes; k++) nIons[k] = 0;
    for (int k = 0; k < nAval; k++) {
        tAvalanche->GetEntry(k);
        // Loop over electrodes
        ionStartPoints = *ionStartPointsInput;
        ionEndPoints = *ionEndPointsInput;
        for (int j = 0; j< (int)ionStartPoints.size(); j++) {
            float zi1 = ionStartPoints[j];
            float zi2 = ionEndPoints[j];
            for (int i = 1; i<nElectrodes+1; i++) { // drift electrode ignored
                if (zi1 < zElectrodes[i] && (zi1>zElectrodes[i+1])) { // l'ion a démarré juste en-dessous de l'électrode d'intérêt
                    hDyingIons[i-1]->Fill(zi2);
                    nIons[i-1]++;
                }
            }
        }
        ionStartPoints.clear();
        ionEndPoints.clear();
    }
    
    TLegend* lgd = new TLegend(0.2, 0.7, 0.9, 0.9);
    for (int j = nElectrodes-1; j>=0; j--) {
        //for (int j = 0; j<nElectrodes; j++) {
        hDyingIons[j]->Scale(1./nIons[j]);
        //hDyingIons[j]->SetMaximum(1.4);
        hDyingIons[j]->SetMaximum(10);
        hDyingIons[j]->GetXaxis()->SetTitle("z (cm)");
        //hDyingIons[j]->SetMaximum(hDyingIons[nElectrodes-1]->GetMaximum());
        hDyingIons[j]->SetLineColor(j+2);
        hDyingIons[j]->SetLineWidth(j+1);
        if (j==nElectrodes-1) hDyingIons[j]->Draw("hist");
        //if (j==0) hDyingIons[j]->Draw("hist");
        else hDyingIons[j]->Draw("hist same");
        lgd->AddEntry(hDyingIons[j], Form("created below %s", electrodeNames[j].c_str()), "l");
    }
    lgd->Draw("same");
}


int FindElectrodeIndex(int modelNum, double zPos) {
    
    vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    double spaceAbove = 0.0030;
    double spaceBelow = 0.0030;
    for (int i = 0; i < (int)zElectrodes.size(); i++) {
        /*
         double spaceAbove = zElectrodes[i-1]-zElectrodes[i];
         double spaceBelow = zElectrodes[i]-zElectrodes[i+1];
         if (zPos < zElectrodes[i]+spaceAbove*factor && zPos > zElectrodes[i]-spaceBelow*factor) return i;
         */
        if (zPos < zElectrodes[i]+spaceAbove && zPos > zElectrodes[i]-spaceBelow) return i;
    }
    return -1;
}


void DrawDyingIonsOverlaid(const int modelNum = 15, TString fSignalName=""){
    
    gStyle->SetLabelSize(.065, "X");
    gStyle->SetLabelOffset(0.01, "X");
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    /*
     //Load geometry parameters
     double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
     LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
     */
    map <string, int, NoSorting> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    int nAval = tAvalanche->GetEntries();
    
    int nElectrodes = electrodeMap.size()-1;    // on s'en fiche des pads
    cout << "There are " << nElectrodes << " electrodes to look at" << endl;
    TH1F* hDyingIons[nElectrodes];
    
    // Initialisation des histos
    int i = 0;
    string electrodeNames[nElectrodes];
    char* elNameChar[nElectrodes];
    map<string, int>::iterator it = electrodeMap.begin();
    //it++;    // ignore drift electrode
    while (it != electrodeMap.end()) {
        if (it->first != "pad") {
            electrodeNames[i] = it->first;
            hDyingIons[i] = new TH1F(Form("hDyingIonsFromBelow%s", electrodeNames[i].c_str()), "Electrode where ions stop", nElectrodes, 0, nElectrodes);
            elNameChar[i] = Form("%s", electrodeNames[i].c_str());
            //strcpy(elNameChar, electrodeNames[i].c_str());
            i++;
        }
        it++;
    }
    
    for (int i = 0; i < nElectrodes; i++) {
        for (int j = 0; j < electrodeNames[i].length(); j++) {
            cout << elNameChar[i][j];
        }
        cout << endl;
    }
    //return;
    
    vector<float> *ionStartPointsInput = 0, *ionEndPointsInput = 0;
    tAvalanche->SetBranchAddress("ionStartPoints", &ionStartPointsInput);
    tAvalanche->SetBranchAddress("ionEndPoints", &ionEndPointsInput);
    
    vector<float> ionStartPoints = {}, ionEndPoints = {};
    int nIons[nElectrodes];
    int nIonsStop[nElectrodes][nElectrodes];
    for (int k = 0; k < nElectrodes; k++) {
        nIons[k] = 0;
        for (int l = 0; l < nElectrodes; l++) {
            nIonsStop[k][l] = 0;
        }
    }
    for (int k = 0; k < nAval; k++) {
        tAvalanche->GetEntry(k);
        // Loop over electrodes
        ionStartPoints = *ionStartPointsInput;
        ionEndPoints = *ionEndPointsInput;
        for (int j = 0; j< (int)ionStartPoints.size(); j++) {
            float zi1 = ionStartPoints[j];
            float zi2 = ionEndPoints[j];
            for (int i = 1; i<nElectrodes+1; i++) { // drift electrode ignored
                if (zi1 < zElectrodes[i] && (zi1>zElectrodes[i+1])) { // l'ion a démarré juste en-dessous de l'électrode d'intérêt
                    //hDyingIons[i-1]->Fill(electrodeNames[i]);
                    //hDyingIons[i-1]->Fill(elNameChar);
                    //hDyingIons[i-1]->Fill(zi2);
                    nIons[i]++;
                    int indexStop = FindElectrodeIndex(modelNum, zi2);
                    if (indexStop < 0) {
                        cout << "Index < 0, problem when drawing transparency" << endl;
                        cout << "Corresponding z is " << zi2 << endl;
                        //return;
                    }
                    nIonsStop[i][indexStop]++;
                }
            }
        }
        ionStartPoints.clear();
        ionEndPoints.clear();
    }
    
    // Combine GEM up and GEM down if there is a GEM
    string sGem = "GEM";
    int newDim = nElectrodes;
    bool gemPresent = false;
    for (int k = 0; k < nElectrodes; k++) {
        if (electrodeNames[k].find(sGem) != string::npos) {
            gemPresent = true;
            newDim = nElectrodes-1;
            break;
        }
    }
    
    
    // Change the dimension of elNameChar and nIonsStop to sum both GEM contributions
    char* elNameCharCopy[newDim];
    string electrodeNamesCopy[newDim];
    int nIonsStopCopy[newDim][newDim];
    // initialisation of nIonsStopCopy et des histos
    for (int k = 0; k < newDim; k++) {
        for (int l = 0; l < newDim; l++) {nIonsStopCopy[k][l]=0;}
    }
    
    // Loop to rename electrodes
    int idx1 = 0;
    for (int k = 0; k < nElectrodes; k++) {
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
    
    // Loop to re-count ions
    idx1 = 0;
    for (int k = 0; k < nElectrodes; k++) {
        int idx2 = 0;
        if (k > 0 && electrodeNames[k].find(sGem)!= std::string::npos && electrodeNames[k-1].find(sGem) != std::string::npos) {idx1--;}
        for (int l = 0; l < nElectrodes; l++) {
            if (l > 0 && electrodeNames[l].find(sGem) != std::string::npos && electrodeNames[l-1].find(sGem) != std::string::npos) {idx2--;}
            nIonsStopCopy[idx1][idx2] += nIonsStop[k][l];
            //cout << "Adding for hist (" << electrodeNamesCopy[idx1] << "," << electrodeNamesCopy[idx2] << ") the values of hist (" << electrodeNames[k] << "," << electrodeNames[l] << ")" << endl;
            idx2++;
        }
        idx1++;
    }
    
    
    // Copy histograms to sum those for the GEM
    TH1F* hDyingIonsCopy[newDim];
    for (int k = 0; k < newDim; k++) {
        hDyingIonsCopy[k] = new TH1F(Form("h2DyingIonsFromBelow%s", electrodeNamesCopy[k].c_str()), "Electrode where ions stop", newDim, 0, newDim);
    }
    
    //return;
    
    int nIonsStopSummed[newDim][newDim];
    double yMin = 100000, yMax = 1;
    for (int i = 0; i<newDim; i++) {
        for (int j = newDim-1; j>=0; j--) {
            nIonsStopSummed[i][j] = 0;
            for (int l = 1; l < i+1; l++) { // ignore the drift
                nIonsStopSummed[i][j] += nIonsStopCopy[l][j];
            }
            //hDyingIons[i]->Fill(elNameChar[j], nIonsStop[i][j]);
            //hDyingIonsCopy[i]->Fill(elNameCharCopy[j], (double)nIonsStopCopy[i][j]/nAval);
            hDyingIonsCopy[i]->Fill(elNameCharCopy[j], (double)nIonsStopSummed[i][j]/nAval);
        }
        if (hDyingIonsCopy[i]->GetMinimum() < yMin) yMin = hDyingIonsCopy[i]->GetMinimum();
        if (hDyingIonsCopy[i]->GetMaximum() > yMax) yMax = hDyingIonsCopy[i]->GetMaximum();
    }
    if (yMin < 0.000001) yMin = 0.001;
    double scale = (yMax-yMin)/10;
    cout << "yMax = " << yMax << endl;
    cout << "scale = " << scale << endl;
    //return;
    scale = 100;
    yMax *= scale;
    yMin = 0.1;
    
    TLegend* lgd = new TLegend(0.2, 0.7, 0.9, 0.9);
    for (int j = newDim-1; j>0; j--) {
        //for (int j = 1; j < newDim; j++) {
        //hDyingIons[j]->Scale(1./nIons[j]);
        //hDyingIons[j]->SetMaximum(10);
        //hDyingIons[j]->GetXaxis()->SetTitle("");
        //hDyingIons[j]->SetMaximum(hDyingIons[nElectrodes-1]->GetMaximum());
        hDyingIonsCopy[j]->SetMinimum(yMin);
        hDyingIonsCopy[j]->SetMaximum(yMax);
        hDyingIonsCopy[j]->SetFillColor(j+1);
        hDyingIonsCopy[j]->SetLineColor(kBlack);
        //hDyingIonsCopy[j]->SetBarWidth(0.98);
        //hDyingIonsCopy[j]->SetBarOffset(0.01);
        //hDyingIonsCopy[j]->SetLineColorAlpha(kWhite, 0.0);
        //hDyingIonsCopy[j]->SetLineWidth(j+1);
        if (j==newDim-1) hDyingIonsCopy[j]->Draw("hist");
        //if (j==1) hDyingIonsCopy[j]->Draw("hist L");
        else hDyingIonsCopy[j]->Draw("hist same");
        if (electrodeNamesCopy[j].find(sGem) != string::npos) lgd->AddEntry(hDyingIonsCopy[j], Form("created in %s", electrodeNamesCopy[j].c_str()), "F");
        else lgd->AddEntry(hDyingIonsCopy[j], Form("created below %s", electrodeNamesCopy[j].c_str()), "F");
    }
    lgd->Draw("same");

    // Draw numbers
    TLatex* txtNions[newDim][newDim];
    for (int i = 1; i < newDim; i++) {  // drift ignored
        for (int j = 0; j < newDim; j++) {
            if ( (double)nIonsStopCopy[i][j]/nAval < 0.1) continue;
            double upperLimit = (double)nIonsStopSummed[i][j]/nAval;
            double lowerLimit = 0.1;
            if (i > 1) lowerLimit = (double)nIonsStopSummed[i-1][j]/nAval;
            if (lowerLimit < 0.1) lowerLimit = 0.1;
            //cout << electrodeNamesCopy[i] << " " << electrodeNamesCopy[j] << " lower limit:" << lowerLimit << " upper limit: " << upperLimit << endl;
            double xPos = newDim-1-j+0.4;
            if ((double)nIonsStopCopy[i][j]/nAval > 1000) xPos = newDim-1-j+0.25;
            else if ((double)nIonsStopCopy[i][j]/nAval > 100)  xPos = newDim-1-j+0.3;
            else if ((double)nIonsStopCopy[i][j]/nAval > 10)  xPos = newDim-1-j+0.35;
            double yPos = TMath::Power(10,(TMath::Log10(lowerLimit)+(TMath::Log10(upperLimit)-TMath::Log10(lowerLimit))*0.45));
            txtNions[i][j] = new TLatex(xPos, yPos, Form("%.0f", (double)nIonsStopCopy[i][j]/nAval));
            if ((double)nIonsStopCopy[i][j]/nAval > 0.1 && (double)nIonsStopCopy[i][j]/nAval < 2) txtNions[i][j] = new TLatex(newDim-1-j+0.3, yPos, Form("%.1f", (double)nIonsStopCopy[i][j]/nAval));
            txtNions[i][j]->Draw();
        }
    }
    
    // Draw lines to seperate the bins
    DrawLines(hDyingIonsCopy[newDim-1]);
}

void DrawDyingIons2d(const int modelNum = 15, TString fSignalName=""){
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    map <string, int, NoSorting> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    int nAval = tAvalanche->GetEntries();
    
    int nElectrodes = electrodeMap.size()-2;    // on s'en fiche de la drift et des pads
    cout << "There are " << nElectrodes << " electrodes to look at" << endl;
    TH2F* hDyingIons[nElectrodes];
    
    // Initialisation des histos de tranparence
    int i = 0;
    string electrodeNames[nElectrodes];
    map<string, int>::iterator it = electrodeMap.begin();
    it++;    // ignore drift electrode
    while (it != electrodeMap.end()) {
        if (it->first != "pad") {
            electrodeNames[i] = it->first;
            hDyingIons[i] = new TH2F(Form("hDyingIonsFromBelow%s", electrodeNames[i].c_str()), "Dying ions = f(Creation of ions)", 100, 0, damp*1.5, 100, 0, ddrift*1.5);
            i++;
        }
        it++;
    }
    
    vector<float> *ionStartPointsInput = 0, *ionEndPointsInput = 0;
    tAvalanche->SetBranchAddress("ionStartPoints", &ionStartPointsInput);
    tAvalanche->SetBranchAddress("ionEndPoints", &ionEndPointsInput);
    
    vector<float> ionStartPoints = {}, ionEndPoints = {};
    
    int nIons[nElectrodes];
    for (int k = 0; k < nElectrodes; k++) nIons[k] = 0;
    for (int k = 0; k < nAval; k++) {
        tAvalanche->GetEntry(k);
        // Loop over electrodes
        ionStartPoints = *ionStartPointsInput;
        ionEndPoints = *ionEndPointsInput;
        for (int j = 0; j< (int)ionStartPoints.size(); j++) {
            float zi1 = ionStartPoints[j];
            float zi2 = ionEndPoints[j];
            for (int i = 1; i<nElectrodes+1; i++) { // drift electrode ignored
                if (zi1 < zElectrodes[i] && (zi1>zElectrodes[i+1])) { // l'ion a démarré juste en-dessous de l'électrode d'intérêt
                    hDyingIons[i-1]->Fill(zi1, zi2);
                    nIons[i-1]++;
                }
            }
        }
        ionStartPoints.clear();
        ionEndPoints.clear();
    }
    
    TLegend* lgd = new TLegend(0.2, 0.7, 0.9, 0.9);
    for (int j = nElectrodes-1; j>=0; j--) {
        hDyingIons[j]->GetXaxis()->SetTitle("z_{Start} (cm)");
        hDyingIons[j]->GetYaxis()->SetTitle("z_{End} (cm)");
        hDyingIons[j]->SetMarkerColor(j+2);
        hDyingIons[j]->SetFillColor(j+2);
        if (j==nElectrodes-1) hDyingIons[j]->Draw("");
        else hDyingIons[j]->Draw("same");
        lgd->AddEntry(hDyingIons[j], Form("created below %s", electrodeNames[j].c_str()), "f");
    }
    lgd->Draw("same");
}

void DrawDyingIons3d(const int modelNum = 15, TString fSignalName=""){
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    
    /*
     gStyle->SetTitleFontSize(.06);
     gStyle->SetTitleSize(.06);
     
     gStyle->SetOptStat(0);
     gStyle->SetTitleFontSize(.05);
     gStyle->SetTitleXSize(.05);
     gStyle->SetTitleYSize(.05);
     gStyle->SetLabelSize(.05, "XY");
     gStyle->SetMarkerStyle(20);
     gStyle->SetMarkerSize(0.3);
     gStyle->SetTextSize(0.05);
     
     gPad->SetLeftMargin(0.15);
     gPad->SetBottomMargin(0.15);
     */
    
    const int unitConversion = 10;    // units go from cm to mm
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    map <string, int, NoSorting> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    int nAval = tAvalanche->GetEntries();
    
    int nElectrodes = electrodeMap.size()-2;    // on s'en fiche de la drift et des pads
    cout << "There are " << nElectrodes << " electrodes to look at" << endl;
    TH3F* hDyingIons[nElectrodes];
    
    // Initialisation des histos de tranparence
    int i = 0;
    string electrodeNames[nElectrodes];
    map<string, int>::iterator it = electrodeMap.begin();
    it++;    // ignore drift electrode
    while (it != electrodeMap.end()) {
        if (it->first != "pad") {
            electrodeNames[i] = it->first;
            hDyingIons[i] = new TH3F(Form("hDyingIonsFromBelow%s", electrodeNames[i].c_str()), "Dying ions", 100, unitConversion*(width/2-pitch), unitConversion*(width/2+pitch), 100, unitConversion*(width/2-pitch), unitConversion*(width/2+pitch), 100, 0, unitConversion*(ddrift*1.1));
            i++;
        }
        it++;
    }
    
    vector<float> *ionStartPointsInput = 0, *ionEndPointsInput = 0, *ionEndPointsInputX = 0, *ionEndPointsInputY = 0;
    tAvalanche->SetBranchAddress("ionStartPoints", &ionStartPointsInput);
    tAvalanche->SetBranchAddress("ionEndPoints", &ionEndPointsInput);
    tAvalanche->SetBranchAddress("ionEndPointsX", &ionEndPointsInputX);
    tAvalanche->SetBranchAddress("ionEndPointsY", &ionEndPointsInputY);
    
    vector<float> ionStartPoints = {}, ionEndPoints = {}, ionEndPointsX = {}, ionEndPointsY = {};
    
    int nIons[nElectrodes];
    for (int k = 0; k < nElectrodes; k++) nIons[k] = 0;
    for (int k = 0; k < nAval; k++) {
        tAvalanche->GetEntry(k);
        // Loop over electrodes
        ionStartPoints = *ionStartPointsInput;
        ionEndPoints = *ionEndPointsInput;
        ionEndPointsX = *ionEndPointsInputX;
        ionEndPointsY = *ionEndPointsInputY;
        for (int j = 0; j< (int)ionStartPoints.size(); j++) {
            float zi1 = ionStartPoints[j];
            float zi2 = ionEndPoints[j];
            float xi2 = ionEndPointsX[j];
            float yi2 = ionEndPointsY[j];
            for (int i = 1; i<nElectrodes+1; i++) { // drift electrode ignored
                if (zi1 < zElectrodes[i] && (zi1>zElectrodes[i+1])) { // l'ion a démarré juste en-dessous de l'électrode d'intérêt
                    hDyingIons[i-1]->Fill(xi2*unitConversion, yi2*unitConversion, zi2*unitConversion);
                    nIons[i-1]++;
                }
            }
        }
        ionStartPoints.clear();
        ionEndPoints.clear();
        ionEndPointsX.clear();
        ionEndPointsY.clear();
    }
    
    TLegend* lgd = new TLegend(0.2, 0.77, 0.9, 0.9);
    for (int j = nElectrodes-1; j>=0; j--) {
        hDyingIons[j]->GetXaxis()->SetTitle("x (mm)");
        hDyingIons[j]->GetYaxis()->SetTitle("y (mm)");
        hDyingIons[j]->GetZaxis()->SetTitle("z (mm)");
        hDyingIons[j]->SetMarkerColor(j+2);
        hDyingIons[j]->SetFillColor(j+2);
        if (j==nElectrodes-1) hDyingIons[j]->Draw("");
        else hDyingIons[j]->Draw("same");
        lgd->AddEntry(hDyingIons[j], Form("created below %s", electrodeNames[j].c_str()), "fP");
    }
    lgd->Draw("same");
}

void DrawNionsInDriftRegion(TString fSignalName=""){
    
    /*
     gStyle->SetTitleFontSize(.06);
     gStyle->SetTitleSize(.06);
     
     gStyle->SetOptStat(0);
     gStyle->SetTitleFontSize(.05);
     gStyle->SetTitleXSize(.05);
     gStyle->SetTitleYSize(.05);
     gStyle->SetLabelSize(.05, "XY");
     gStyle->SetMarkerSize(0.3);
     gStyle->SetMarkerStyle(20);
     gStyle->SetTextSize(0.05);
     */
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    int ionBackNum;
    //int amplificationElectrons;
    vector <Int_t> *nWinnersVec = nullptr;
    vector <Int_t> gainVec = {};
    tAvalanche->SetBranchAddress("ionBackNum", &ionBackNum);
    //tAvalanche->SetBranchAddress("amplificationElectrons", &amplificationElectrons);
    tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVec);
    int nMax = tAvalanche->GetMaximum("ionBackNum");
    
    
    TH1F* hIonsNumber = new TH1F("hIonsNumber", "Number of ions in the drift region", nMax, 0, nMax);
    int nAval = tAvalanche->GetEntries();
    
    for (int k = 0; k < nAval; k++) {
        tAvalanche->GetEntry(k);
        gainVec = *nWinnersVec;
        //if (amplificationElectrons>1) hIonsNumber->Fill(ionBackNum);
        if (gainVec[0]>1) hIonsNumber->Fill(ionBackNum);
    }
    hIonsNumber->GetXaxis()->SetTitle("Number of ions");
    hIonsNumber->Draw("");
    
    
}

void DrawTransparency(const int modelNum = 15, TString fSignalName="") {
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    
    Int_t elMax = GetMaxAmp(*tAvalanche);    // must be called before setting the addresses to the tree
    
    map <string, int, NoSorting> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    /*
     for (int i = 0; i<zElectrodes.size(); i++) {cout << zElectrodes[i] << endl;}
     return;
     */
    
    int nAval = tAvalanche->GetEntries();
    
    // First draw total transparency
    //Int_t nAmplification;
    vector <Int_t> *nWinnersVec = nullptr;
    vector <Int_t> gainVec = {};
    int winners = 0;
    //tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
    tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVec);
    TH1F* hTransparencyTotal = new TH1F("hTransparency", "Total transparency", 4, -1, 3);
    for (int k = 0; k<nAval; k++) {
        tAvalanche->GetEntry(k);
        gainVec = *nWinnersVec;
        if (gainVec[0] > 1) {
            //if (nAmplification>1) {
            hTransparencyTotal->Fill(1);
            winners++;
        }
        else {hTransparencyTotal->Fill(0);}
    }
    hTransparencyTotal->Scale(1./nAval);
    hTransparencyTotal->SetMaximum(1.2);
    // This is just a boolean histogram, don't draw it in the end
    //hTransparencyTotal->Draw("hist");
    
    // Write text with the value of transparency
    Double_t totalTransp = (double)winners/nAval;
    
    
    // Now draw transparency of all electrodes
    int nElectrodes = electrodeMap.size()-2;    // on s'en fiche de la drift et des pads
    //cout << "There are " << nElectrodes << " electrodes to look at" << endl;
    
    // Initialisation des histos de tranparence et de gain
    int i = 0;
    string electrodeNames[nElectrodes];
    TH1F* hTransparency[nElectrodes];
    TH1F* hGains[nElectrodes];
    
    map<string, int>::iterator it = electrodeMap.begin();
    it++;    // ignore drift electrode
    while (it != electrodeMap.end()) {
        if (it->first != "pad") {
            electrodeNames[i] = it->first;
            hTransparency[i] = new TH1F(Form("transparency%s", electrodeNames[i].c_str()), Form("Transparency of %s", electrodeNames[i].c_str()), 4, -1, 3);
            hGains[i] = new TH1F(Form("gain%s", electrodeNames[i].c_str()), Form("Gain of %s", electrodeNames[i].c_str()), elMax, 0, 2*elMax);
            i++;
        }
        it++;
    }
    
    //return;
    
    vector<float> *electronStartPointsInput = 0, *electronEndPointsInput = 0;
    //vector<float> *ionStartPointsInput = 0, *ionEndPointsInput = 0;
    
    tAvalanche->SetBranchAddress("electronStartPoints", &electronStartPointsInput);
    tAvalanche->SetBranchAddress("electronEndPoints", &electronEndPointsInput);
    //tAvalanche->SetBranchAddress("ionStartPoints", &ionStartPointsInput);
    //tAvalanche->SetBranchAddress("ionEndPoints", &ionEndPointsInput);
    
    vector<float> electronStartPoints = {}, electronEndPoints = {};
    vector<float> ionStartPoints = {}, ionEndPoints = {};
    
    int nElectrons[nElectrodes];        // number of electrons that exist above the electrode of interest
    double nElCreated[nElectrodes];        // average number of electrons created below the electrode of interest
    // initialization
    for (int i = 0; i<nElectrodes; i++) {
        nElectrons[i] = 0;
        nElCreated[i] = 0;
    }
    //int nElectronsTotal = 0;
    for (int k = 0; k < nAval; k++) {
        tAvalanche->GetEntry(k);
        
        //nElectronsTotal += electronStartPoints.size();
        // Loop over electrodes
        electronStartPoints = *electronStartPointsInput;
        electronEndPoints = *electronEndPointsInput;
        //ionStartPoints = *ionStartPointsInput;
        //ionEndPoints = *ionEndPointsInput;
        for (int j = 0; j< (int)electronStartPoints.size(); j++) {
            float ze1 = electronStartPoints[j];
            float ze2 = electronEndPoints[j];
            for (int i = 1; i < nElectrodes+1; i++) { // drift electrode ignored
                // First, transparency
                if (ze1 > zElectrodes[i] && ze2 < (zElectrodes[i-1]- (zElectrodes[i-1]-zElectrodes[i])/3) ) { // il y a bien un électron au-dessus de l'électrode d'intérêt et il n'est pas mort au stade encore au-dessus
                    nElectrons[i-1]++;
                    if (ze2 < zElectrodes[i] - (zElectrodes[i]-zElectrodes[i+1])/2) hTransparency[i-1]->Fill(1);
                    else hTransparency[i-1]->Fill(0);
                }
                // Now gain
                if (ze1 < zElectrodes[i] && ze1 > zElectrodes[i+1] ) { // un électron a été créé en-dessous de l'électrode d'intérêt et au-dessus de l'électrode de dessous
                    nElCreated[i-1]++;
                    //cout << "ok" << endl;
                }
            }
        }
        electronStartPoints.clear();
        electronEndPoints.clear();
        //ionStartPoints.clear();
        //ionEndPoints.clear();
    }
    // normalize with the number of entries
    for (int i = 0; i < nElectrodes; i++) {nElCreated[i] /= nAval;}
        
    
    //cout << "nElectronsTotal = " << nElectronsTotal << endl;
    double transp[nElectrodes];
    for (int j = 0; j<nElectrodes; j++) {
        int bin = hTransparency[j]->GetXaxis()->FindBin(1);
        cout << "bin = " << bin << endl;
        cout << "nElectrons[j] = " << nElectrons[j] << endl;
        if (nElectrons[j] == 0) {cout << "\n\n\nnElectrons = 0!!\n\n\n" << endl; return;}
        transp[j] = (double)hTransparency[j]->GetBinContent(bin)/nElectrons[j];
        //cout << "nElectrons[j] between ze1>" << zElectrodes[j+1] << " and ze2 < " << zElectrodes[j] << " = " << nElectrons[j] << endl;
        //hTransparency[j]->Scale(1./nElectronsTotal);
        hTransparency[j]->Scale(1./nElectrons[j]);
        hTransparency[j]->SetMaximum(1.2);
        //hTransparency[j]->Draw("hist same");
        TLatex* txt = new TLatex(0.2,0.9-(j+1)*0.07,Form("#bf{%s transparency = %.1f %s}", electrodeNames[j].c_str(), transp[j]*100, "%"));
        txt->Draw("same");
    }
    
    //TLatex* txttr = new TLatex(-0.5,0.9*hTransparencyTotal->GetMaximum(),Form("#bf{Full detector transparency = %.1f %s}", totalTransp*100, "%"));
    TLatex* txttr = new TLatex(0.2,0.9,Form("#bf{Full detector transparency = %.1f %s}", totalTransp*100, "%"));
    txttr->Draw("same");
    
    // Draw gain texts
    double nElAbove = 1;    // one electron above the drift electrode at first
    //double nElCreated[nElectrodes];
    double gain[nElectrodes];
    double yMax = 0.8-0.07*(nElectrodes+1);
    for (int j = 0; j<nElectrodes; j++) {
        //double gain = (double)hGains[j]->Integral()/hGains[j]->GetEntries();
        /*
        nElCreated[j] = 0;
        for (int l = 0; l < hGains[j]->GetNbinsX(); l++) {
            if (hGains[j]->GetBinContent(l) > 0) {
                double val = hGains[j]->GetXaxis()->GetBinCenter(l);
                nElCreated[j] += hGains[j]->GetBinContent(l)*val;
            }
        }
        nElCreated[j] /= hGains[j]->GetEntries();   // on normalise avec le nombre d'événements
         */
        //cout << endl << "index " << j << " corresponding to " << electrodeNames[j] << endl;
        //cout << "number of electrons created = " << nElCreated[j] << endl;
        if (j>0) {nElAbove *= gain[j-1];}
        //cout << "number of electrons above = " << nElAbove << endl;
        nElAbove *= transp[j];
        gain[j] = nElCreated[j] / nElAbove;
        //cout << "gain = nElCreated[j] / (nElAbove*transparency) = " << gain[j] << endl;
        double gainError = 1/TMath::Sqrt(hGains[j]->GetEntries());
        //TLatex* txt = new TLatex(0.2,yMax-j*0.07,Form("#bf{Gain of %s = %.1f #pm %.1f}", electrodeNames[j].c_str(), gain[j], gainError));
        TLatex* txt = new TLatex(0.2,yMax-j*0.07,Form("#bf{Gain of %s = %.1f}", electrodeNames[j].c_str(), gain[j]));
        txt->Draw("same");
    }
    
    /*
    TCanvas* cvTest = new TCanvas();
    cvTest->Divide(5,5);
    for (int j = 0; j<nElectrodes; j++) {
        cvTest->cd(j+1);
        hTransparency[j]->Draw();
    }
    TCanvas* cvTest2 = new TCanvas();
    cvTest2->Divide(5,5);
    for (int j = 0; j<nElectrodes; j++) {
        cvTest2->cd(j+1);
        hGains[j]->Draw();
    }
    cvTest->SaveAs("Figures/test-transp.pdf");
    cvTest2->SaveAs("Figures/test-gains.pdf");
     */
}


int Transparency() {
    /* Main function, for testing */
    
    //______________________
    // variables
    string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 10;
    vector <Int_t> hvList = {330, 410, 530, 730, 850};
    //____________________
    time_t t0 = time(NULL);
    
    gStyle->SetTitleFontSize(.06);
    gStyle->SetTitleSize(.06);
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    gStyle->SetMarkerSize(0.3);
    gStyle->SetTextSize(0.05);
    gStyle->SetLabelOffset(0.01);
    
    
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    TString fSignalName = path + "signal";
    for (int k = 0; k < (int)hvList.size(); k++) {fSignalName += Form("-%d", hvList[k]);}
    fSignalName += ".root";
    cout << fSignalName << endl;
    
    
    TString outputName = Form("Figures/model%d/transparency-2d-%s", modelNum, gasName.c_str());
    for (int k = 0; k< (int)hvList.size(); k++) {outputName += Form("-%d", hvList[k]);}
    outputName+=".pdf";
    
    
    TCanvas* cv = new TCanvas("cv","cv", 1000, 800);
    //DrawTransparency(modelNum, fSignalName);
    DrawDyingIons2d(modelNum, fSignalName);
    cv->SaveAs(outputName);
    
    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    return 0;
    
}
