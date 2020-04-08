#include "TNtuple.h"
#include "TString.h"
#include "TDatime.h"
#include "TFile.h"
#include "TColor.h"
#include "TStyle.h"
#include "TColor.h"
#include <iostream>
#include <fstream>
#include <dirent.h>

#include "_Utils.C"

// to draw i(t) for a series of files on the 4 electrodes + the histograms of i_mon

// calcul d'ibf à faire ici
/*
// A faire :
 - récupérer les courants de la drift et du pad (avec la condition que le courant des ions n'est pas nul, sinon on met dans des histogrammes des valeurs qui ne correspondent pas à des événements)
 - les mettre dans des histogrammes
 - fitter ces histogrammes pour en tirer l'ibf
 
 --> extraire aussi la valeur de courant moyenne
 --> évaluer la proba de pile up
 
 */

double GetUpperLimit(double x) {
    if (x > 0) return x*1.1;
    else return x/1.1;
}

double GetLowerLimit(double x) {
    if (x > 0) return x/1.1;
    else return x*1.1;
}

Double_t FitGauss( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
return  par[2]*TMath::Gaus( x[0], par[0], par[1]); }

TF1* FitCurrent(TH1F* hCurrent, Double_t fitWidth) {
    Int_t iBinMax = hCurrent->GetMaximumBin();
    Double_t xMax = hCurrent->GetXaxis()->GetBinCenter( iBinMax );
    
    std::cout << "xMax = " << xMax << std::endl;
    std::cout << "maximum = " << hCurrent->GetMaximum() << std::endl;
    
    //Int_t fitRangeMin = xMax - 0.6 * hCurrent->GetRMS();
    //Int_t fitRangeMax = xMax + 0.6 * hCurrent->GetRMS();
    Int_t fitRangeMin = xMax - fitWidth;
    Int_t fitRangeMax = xMax + fitWidth;

    TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
    f->SetParNames("Mean", "Sigma", "Amplitude");
    f->SetParameters(xMax, hCurrent->GetRMS(), hCurrent->GetMaximum());
    
    hCurrent->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
    return f;
}


void GetIbf() {
    

    gStyle->SetTitleFontSize(.06);
    gStyle->SetTitleSize(.06);

    /*
     gStyle->SetOptStat(0);
     gStyle->SetTitleFontSize(.05);
     gStyle->SetTitleXSize(.05);
     gStyle->SetTitleYSize(.05);
     gStyle->SetLabelSize(.05, "XY");
     */
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    //____________________

    time_t t0 = time(NULL);
    gStyle->SetOptStat(0);
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);

     
    // Get number of files to look at
    Int_t num = 1;
    //Int_t num = GetNumberOfFiles(path, "signal");
    /*Int_t num2 = int(num/2.);
    if (num/2.> num2) num2+=1;
    TCanvas* c2 = new TCanvas("c2");
    c2->Divide(2, num2);
    */
    Int_t electrodeNum = GetElectrodeNum(modelNum);
    
    Double_t ibfList[num], ibfErrorList[num];
    
    
    /*
    for (unsigned int k = 0; k < num; ++k) {
        Int_t hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvGemDown = 0, hvGemUp = 0, hvDrift = 0;
        TString signalFileName, fOutputName;
        if (modelNum == 1) {
            hvMesh = 340+20*k;
            hvDrift = 540+20*k;
            signalFileName = path + Form("signal-%d-%d.root", hvMesh, hvDrift);
            fOutputName = path + Form("amp-%d-%d.root", hvMesh, hvDrift);
            //
            signalFileName = path + "signal-380-580.root";
            fOutputName = path + "amp-380-580.root";
        }
        else if (modelNum >= 2 && modelNum < 5) {
            hvDmDown = 300;
            hvDmUp = 600;
            hvDrift = 800;
            signalFileName = path + Form("signal-%d-%d-%d.root", hvDmDown, hvDmUp, hvDrift);
            fOutputName = path + Form("amp-%d-%d-%d.root", hvDmDown, hvDmUp, hvDrift);
        }
        else if (modelNum >= 5 && modelNum < 8) {
            hvMesh = 300;
            hvDmDown = 400;
            hvDmUp = 700;
            hvDrift = 900;
            signalFileName = path + Form("signal-%d-%d-%d-%d.root", hvMesh, hvDmDown, hvDmUp, hvDrift);
            fOutputName = path + Form("amp-%d-%d-%d-%d.root", hvMesh, hvDmDown, hvDmUp, hvDrift);
        }
        else if (modelNum >= 8 && modelNum < 14) {
            hvMesh = 350;
            hvDmDown = 430;
            hvDmUp = 630;
            hvDrift = 750;
            signalFileName = path + Form("signal-%d-%d-%d-%d.root", hvMesh, hvGemDown, hvGemUp, hvDrift);
            fOutputName = path + Form("amp-%d-%d-%d-%d.root", hvMesh, hvGemDown, hvGemUp, hvDrift);
        }
     */
        // import data
        //TFile* fSignal = TFile::Open(signalFileName, "READ");
        //TFile* fAmp = new TFile(fOutputName, "RECREATE");
    
    //TFile* fSignal = TFile::Open("rootFiles/Ar-iC4H10/model1/feSignal-340-540.root", "READ");
    TFile* fSignal = TFile::Open("rootFiles/Ar-iC4H10/model1/signal-340-540.root", "READ");

    
    std::map <std::string, int> electrodeMap;
    LoadElectrodeMap(modelNum, electrodeMap);

    Double_t fctPad = 0.;
    //Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
    TTree* tSignalPad = (TTree*)fSignal->Get(Form("tSignal_%d", electrodeMap["pad"]));
    //tSignal->SetBranchAddress("time", &ft);
    tSignalPad->SetBranchAddress("totalCurrent", &fctPad);
    //tSignal->SetBranchAddress("electronCurrent", &fce);
    //tSignal->SetBranchAddress("ionCurrent", &fci);
    
    Double_t fctDrift = 0.;
    TTree* tSignalDrift = (TTree*)fSignal->Get(Form("tSignal_%d", electrodeMap["drift"]));
    tSignalDrift->SetBranchAddress("totalCurrent", &fctDrift);
    
    
    int nEntries = tSignalPad->GetEntries();
    if (tSignalPad->GetEntries() != tSignalDrift->GetEntries() ) {
        std::cout << "not the same number of entries in tree pad and tree drift" << std::endl;
        return;
    }
    
    double iDriftMax = tSignalDrift->GetMaximum("totalCurrent");
    double iDriftMin = tSignalDrift->GetMinimum("totalCurrent");
    double iPadMax = tSignalPad->GetMaximum("totalCurrent");
    double iPadMin = tSignalPad->GetMinimum("totalCurrent");
    
    //TH1F* hCurrentPad = new TH1F("hCurrentPad", "Currents in the pad", 400, GetLowerLimit(iPadMin), GetUpperLimit(iPadMax));
    TH1F* hCurrentPad = new TH1F("hCurrentPad", "Currents in the pad", 10000, iPadMin-20, iPadMax+20);
    TH1F* hCurrentDrift = new TH1F("hCurrentDrift", "Currents in the drift electrode", 5000, iDriftMin-1, iDriftMax+1);
    for (int k=0; k<nEntries; k++) {
        tSignalPad->GetEntry(k);
        tSignalDrift->GetEntry(k);
        if (abs(fctPad) > 10) hCurrentPad->Fill(fctPad);
        if (abs(fctDrift) > 0) hCurrentDrift->Fill(fctDrift);
    }
    
    TCanvas* cv = new TCanvas("cv", "cv", 800, 400);
    cv->Divide(2);
    
    cv->cd(1);
    TF1* fPad = FitCurrent(hCurrentPad, 3.);
    hCurrentPad->GetXaxis()->SetRangeUser(fPad->GetParameter(0)-20, fPad->GetParameter(0)+20);
    hCurrentPad->Draw();
    fPad->Draw("same");
    cv->cd(2);
    TF1* fDrift = FitCurrent(hCurrentDrift, 0.01);
    //hCurrentDrift->GetXaxis()->SetRangeUser(fDrift->GetParameter(0)-20, fDrift->GetParameter(0)+20);
    hCurrentDrift->Draw();
    //fDrift->SetLineColor(kGreen);
    fDrift->Draw("same");
    
/*
            TLegend* legend = new TLegend(0.55,0.7,0.89,0.8);
            legend->AddEntry(hImonWith[j],"With source","l");
            legend->AddEntry(hImonWithout[j],"Without source","l");
            legend->Draw();
 */
    cv->SaveAs("Figures/currents-signal.pdf");
    
    ibfList[0] = -fDrift->GetParameter(0)/fPad->GetParameter(0);
    std::cout << "ibf = " << ibfList[0] << "%" << std::endl;
}
