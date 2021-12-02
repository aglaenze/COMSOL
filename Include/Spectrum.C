#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>

#include <TCanvas.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>

#include "Functions.C"
#include "Utils.C"
#include "Geometry.C"
#include "Transparency.C"

using namespace std;

Double_t ln(Double_t x) {
    return TMath::Log(x);
}

void WriteInfo(TF1* f) {
    
    // Write gain, sigma and res = sigma/gain on the plot
    Double_t gain = f->GetParameter(0);
    Double_t gainError = f->GetParError(0);
    Double_t sigma = abs(f->GetParameter(1));
    Double_t sigmaError = f->GetParError(1);
    Double_t resolution = f->GetParameter(1)/f->GetParameter(0);
    Double_t resolutionError = resolution * TMath::Sqrt( Square(f->GetParError(0)/f->GetParameter(0)) + Square(f->GetParError(1)/f->GetParameter(1)) );
    Double_t fwhm = 2*TMath::Sqrt(2* TMath::Log(2))*sigma;
    Double_t fwhmError = 2*TMath::Sqrt(2* TMath::Log(2))*sigmaError;
    /*
     Double_t fwhm = resolution * 2* ln(2);
     Double_t fwhmError = resolutionError * 2* ln(2);
     */
    
    double xPos = f->GetParameter(0)*0.2;
    TLatex* txtGain = new TLatex(xPos, 1.2, Form("#bf{Gain = %.1f #pm %.1f}", gain, gainError));
    TLatex* txtSigma = new TLatex(xPos, 1.1, Form("#bf{#sigma = %.1f #pm %.1f}", sigma, sigmaError));
    
    
    std::string percent = "%";
    TLatex* txtRes = new TLatex(xPos, 1.0, Form("#bf{#sigma/Mean = (%.1f #pm %.1f) %s}", abs(resolution)*100, abs(resolutionError)*100, percent.c_str()));
    TLatex* txtFwhm = new TLatex(xPos, 1.0, Form("#bf{FWHM = %.1f #pm %.1f}", abs(fwhm), abs(fwhmError)));
    
    txtGain->Draw("same");
    txtSigma->Draw("same"); //txtRes->Draw("same");
    txtFwhm->Draw("same");
}

void DrawAmplificationElectrons(string gasName = "Ar-iC4H10", TString fSignalName="", bool useFeSource = false) {
    /* Draw the gain */
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
    Int_t nAvalanche = tAvalanche->GetEntries();
    
    if (useFeSource) {
        cout << "There are " << (int)tAvalanche->GetEntries() << " entries in this file" << endl;
    }
    
    Int_t nPrimaryTh = GetPrimary(gasName);
    
    // Create histogram of amplification electrons
    Int_t elMax = GetMaxAmp(*tAvalanche);    // must be called before setting the addresses to the tree
    if (useFeSource) elMax = int((double)elMax/nPrimaryTh);
    
    /*
     Int_t nAmplification;
     tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
     */
    vector<Int_t> *neAvalVecIn = nullptr, *nWinnersVecIn = nullptr;
    tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVecIn);
    tAvalanche->SetBranchAddress("avalancheSize", &neAvalVecIn);
    
    vector<Int_t> nWinnersVec = {}, neAvalVec = {};
    TH1F* hAmplification = new TH1F("hAmplification", "Number of amplification electrons", elMax, 0, 2*elMax);
    const int nEntries = (int)tAvalanche->GetEntries();
    for (int k = 0; k< (int)tAvalanche->GetEntries(); k++) {
        tAvalanche->GetEntry(k);
        nWinnersVec = *nWinnersVecIn;
        neAvalVec = *neAvalVecIn;
        int nAmplification = 0;
        //cout << nWinnersVec.size() << endl;
        for (int l = 0; l< (int)nWinnersVec.size(); l++) {nAmplification += nWinnersVec[l];}
        //cout << nAmplification << endl;
        if (useFeSource) nAmplification /= (double)nPrimaryTh;
        if (nAmplification > 5) hAmplification->Fill(nAmplification);
        //cout << (double)nAmplification/nPrimaryTh << endl;
        nWinnersVec.clear();
        neAvalVec.clear();
    }
    while (hAmplification->GetMaximum() < 30) hAmplification->Rebin(2);
    if ( nEntries > 1000) {while (hAmplification->GetMaximum() < 70) hAmplification->Rebin(2);}
    if ( nEntries > 5000) {while (hAmplification->GetMaximum() < 120) hAmplification->Rebin(2);}
    if ( nEntries > 10000) {while (hAmplification->GetMaximum() < 200) hAmplification->Rebin(2);}
    hAmplification->Scale(1/hAmplification->GetMaximum());
    hAmplification->SetMaximum(1.3);
    hAmplification->SetMinimum(0);
    hAmplification->SetLineColor(kBlue);
    
    TF1* fAmplification = GetFitCurve(hAmplification, useFeSource);    // if Fe source, fit = gauss; else fit = landau
    fAmplification->SetLineColor(kBlue);
    
    Int_t iBinMax = hAmplification->GetMaximumBin();
    Double_t xMax = hAmplification->GetXaxis()->GetBinCenter( iBinMax );
    hAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hAmplification->GetRMS());
    hAmplification->GetXaxis()->SetTitle("Number of electrons");
    
    if (xMax > 8000) hAmplification->GetXaxis()->SetMaxDigits(3);
    
    hAmplification->Draw("hist");
    fAmplification->Draw("same");
    
    // Add text to frame
    TString txt = Form("Number of electrons --> Gain = %.0f #pm %.3f", fAmplification->GetParameter(0), fAmplification->GetParError(0));
    TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
    legend->AddEntry(fAmplification,txt,"l");
    legend->SetTextSize(0.04);
    //legend->Draw("same");
    
    // Write gain, sigma and res = sigma/gain on the plot
    WriteInfo(fAmplification);
}


void DrawFeConvolution(TString fConvolutedName, Double_t& gain, Double_t& gainError) {
    // First draw convoluted spectrum
    TFile* fConvoluted = TFile::Open(fConvolutedName, "READ");
    
    TH1F* hFeAmplification = (TH1F*)fConvoluted->Get("hFeAmplification"); // total number of electrons
    while (hFeAmplification->GetMaximum() < 100) hFeAmplification->Rebin(2);
    
    hFeAmplification->Scale(1/hFeAmplification->GetMaximum());
    hFeAmplification->SetMaximum(1.3);
    //hFeAmplification->GetXaxis()->SetRangeUser(2, 10000);
    hFeAmplification->SetLineColor(kBlue);
    hFeAmplification->SetTitle("Gain with Fe source");
    hFeAmplification->GetXaxis()->SetTitle("Normalised number of electrons");
    
    TF1* f = GetFitCurve(hFeAmplification);
    f->SetLineColor(kBlue);
    
    Int_t iBinMax = hFeAmplification->GetMaximumBin();
    Double_t xMax = hFeAmplification->GetXaxis()->GetBinCenter( iBinMax );
    hFeAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
    
    if (xMax > 8000) hFeAmplification->GetXaxis()->SetMaxDigits(3);
    
    hFeAmplification->Draw("hist");
    f->Draw("same");
    
    // Write gain, sigma and res = sigma/gain on the plot
    WriteInfo(f);
    
    /*
     TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
     legend->AddEntry(f,txt,"l");
     legend->AddEntry(f2,txt2,"l");
     legend->SetTextSize(0.04);
     legend->Draw("same");
     */
}

void DrawFeConvolution(TString fConvolutedName="") {
    Double_t gain = 0, gainError = 0;
    DrawFeConvolution(fConvolutedName, gain, gainError);
    return;
}

void DrawFeChargeConvolution(int modelNum, TString fConvolutedName, string readout, Double_t& gain, Double_t& gainError) {
    
    std::map <std::string, int, NoSorting> electrode;
    LoadElectrodeMap(modelNum, electrode);
    int readoutElectrode = 0;
    std::map<std::string, int>::iterator it = electrode.begin();
    for (it=electrode.begin(); it!=electrode.end(); ++it) {
        //std::cout << it->first << " => " << it->second << '\n';
        if (it->first == readout) readoutElectrode = it->second;    // could be mesh, it depends on where you want to read
    }
    if (readoutElectrode == 0) {
        std::cout << "Did not find mesh electrode, using pad instead" << std::endl;
        readout = "pad";
        for (it=electrode.begin(); it!=electrode.end(); ++it) {
            if (it->first == readout) readoutElectrode = it->second;    // could be mesh, it depends on where you want to read
        }
        if (readoutElectrode == 0) {
            std::cout << "Did not find pad electrode" << std::endl;
            return;
        }
    }
    
    TFile* fConvoluted = TFile::Open(fConvolutedName, "READ");
    TH1F* hFeCharge = (TH1F*)fConvoluted->Get(Form("hFeCharge_%d", readoutElectrode));
    while (hFeCharge->GetMaximum() < 100) hFeCharge->Rebin(2);
    //hFeAmplification->Rebin(8);
    
    hFeCharge->Scale(1/hFeCharge->GetMaximum());
    hFeCharge->SetMaximum(1.3);
    hFeCharge->SetLineColor(kRed);
    hFeCharge->SetLineColor(kBlue);
    hFeCharge->SetTitle("Gain with Fe source");
    hFeCharge->GetXaxis()->SetTitle(Form("Normalised number of charges on the %s", readout.c_str()));
    
    
    TF1* f = GetFitCurve(hFeCharge);
    f->SetLineColor(kRed);
    
    Int_t iBinMax = hFeCharge->GetMaximumBin();
    Double_t xMax = hFeCharge->GetXaxis()->GetBinCenter( iBinMax );
    hFeCharge->GetXaxis()->SetRangeUser(0, xMax + 3*hFeCharge->GetRMS());
    
    if (xMax > 8000) hFeCharge->GetXaxis()->SetMaxDigits(3);
    
    hFeCharge->Draw("hist same");
    f->SetLineColor(kBlue);
    f->Draw("same");
    TString txt2 = Form("Induced charge --> Gain = %.0f #pm %.3f", f->GetParameter(0), f->GetParError(0));
    /*
     TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
     legend->AddEntry(f,txt,"l");
     legend->AddEntry(f2,txt2,"l");
     legend->SetTextSize(0.04);
     legend->Draw("same");
     */
    
    WriteInfo(f);
}

void DrawFeChargeConvolution(int modelNum = 0, TString fConvolutedName="", string readout = "") {
    Double_t gain = 0, gainError = 0;
    DrawFeChargeConvolution(modelNum, fConvolutedName, readout, gain, gainError);
    return;
}

void DrawGains(int modelNum = 0, TString fConvolutedName="") {
    /* Not used in the end */
    
    TFile* fConvoluted = TFile::Open(fConvolutedName, "READ");
    
    std::map <std::string, int, NoSorting> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);
    std::vector<double> zElectrodes = {};
    LoadParameters(modelNum, zElectrodes);
    
    int nElectrodes = electrodeMap.size()-2;    // on s'en fiche de la drift et des pads
    std::cout << "There are " << nElectrodes << " electrodes to look at" << std::endl;
    
    string electrodeNames[nElectrodes];
    double gains[nElectrodes-1], gainErrors[nElectrodes];
    map<string, int>::iterator it = electrodeMap.begin();
    it++;    // ignore drift electrode
    int k = 0;
    while (it != electrodeMap.end()) {
        if (it->first != "pad") {
            electrodeNames[k] = it->first;
            int elIndex = it->second;
            // Get gain for each electrode
            TH1F* hFeCharge = (TH1F*)fConvoluted->Get(Form("hFeCharge_%d", elIndex));
            while (hFeCharge->GetMaximum() < 100) hFeCharge->Rebin(2);
            TF1* f2 = GetFitCurve(hFeCharge);
            gains[k] = f2->GetParameter(0);
            gainErrors[k] = f2->GetParError(0);
            k++;
        }
        it++;
    }
    
    double yMax = 0.8-0.07*(nElectrodes+1);
    for (int j = 0; j<nElectrodes; j++) {
        TLatex* txt = new TLatex(0.2,yMax-j*0.07,Form("#bf{Gain of %s = %.1f #pm %.1f}", electrodeNames[j].c_str(), gains[j], gainErrors[j]));
        txt->Draw("same");
    }
    
}

int Spectrum(int modelNum, std::string gasName, std::vector<int> hvList, bool useFeSource) {
    
    time_t t0 = time(NULL);
    LoadStyle();
    
    int electrodeNum = GetElectrodeNum(modelNum);
    if ((int)hvList.size() != electrodeNum-1) {cout << "Wrong hv input" << endl; return 0;}
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    // input and output files
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    TString fConvolutedName = path + "fe-spectrum-convoluted";
    TString fSignalName = path;
    if (useFeSource) fSignalName += "fe-signal-noibf";
    else fSignalName += "signal-noibf";
    //TString fSignalName = path + "signal";
    TString outputName = Form("Figures/model%d/spectra-%s", modelNum, gasName.c_str());
    for (int k = 0; k< (int)hvList.size(); k++) {
        fConvolutedName += Form("-%d", hvList[k]);
        fSignalName += Form("-%d", hvList[k]);
        outputName += Form("-%d", hvList[k]);
    }
    fConvolutedName += ".root";
    fSignalName += ".root";
    outputName+=".pdf";
    
    cout << fSignalName << endl;
    if (useFeSource) {cout << endl << "USING FE SOURCE" << endl << endl;}
    
    
    TCanvas* cv = new TCanvas("cv","cv", 1200, 1000);
    cv->Divide(2);
    
    cv->cd(1);
    DrawAmplificationElectrons(gasName, fSignalName, useFeSource);
    
    cv->cd(2);
    DrawFeConvolution(fConvolutedName);
    
    /*
     TLegend* lgd = new TLegend(0.3, 0.8, 0.9, 0.9);
     lgd->AddEntry(hFeAmplification, "Convolution with Fe spectrum", "l");
     //lgd->AddEntry(hAmplification, "Fe source simulated", "l");
     lgd->Draw();
     */
    
    cv->SaveAs(outputName);
    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    return 0;
}


