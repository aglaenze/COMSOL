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

#include "Spectrum.C"

using namespace std;

void DrawIbf(int modelNum = 0, TString fSignalName="") {
    /* Draw the ibf */
    
    std::map <std::string, int, NoSorting> electrode;
    LoadElectrodeMap(modelNum, electrode);
    int readoutElectrode = 0;
    int driftElectrode = 0;
    
    std::map<std::string, int>::iterator it = electrode.begin();
    for (it=electrode.begin(); it!=electrode.end(); ++it) {
        //std::cout << it->first << " => " << it->second << '\n';
        if (it->first == "mesh") readoutElectrode = it->second;    // could be mesh, it depends on where you want to read
        else if (it->first == "drift") driftElectrode = it->second;
    }
    if (readoutElectrode == 0) {
        std::cout << "Did not find mesh electrode, trying pad" << std::endl;
        for (it=electrode.begin(); it!=electrode.end(); ++it) {
            if (it->first == "pad") readoutElectrode = it->second;
        }
    }
        
    if (readoutElectrode == 0 || driftElectrode == 0) {std::cout << "Did not find drift or pad electrode" << std::endl; return;}
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
    Int_t nAvalanche = tAvalanche->GetEntries();
    
    /*
     Int_t nAmplification;
     tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
     */
    vector<Int_t> *neAvalVecIn = nullptr, *nWinnersVecIn = nullptr;
    tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVecIn);
    tAvalanche->SetBranchAddress("avalancheSize", &neAvalVecIn);
    Int_t ni = 0, ionBackNum = 0;
    tAvalanche->SetBranchAddress("ionNum", &ni);
    tAvalanche->SetBranchAddress("ionBackNum", &ionBackNum);
    vector<Int_t> nWinnersVec = {}, neAvalVec = {};
    
    TH1F* hIbf = new TH1F("ibf", "ibf", 10000, 0, 100);
    for (int k = 0; k<nAvalanche; k++) {
        tAvalanche->GetEntry(k);
        //if (nAmplification>1) {hIbf->Fill((double)ionBackNum/ni*100.);}
        nWinnersVec = *nWinnersVecIn;
        neAvalVec = *neAvalVecIn;
        int nAmplification = 0;
        for (int l = 0; l< (int)nWinnersVec.size(); l++) {nAmplification += nWinnersVec[l];}
        if (nAmplification > 2) {hIbf->Fill((double)ionBackNum/ni*100.);}
        nWinnersVec.clear();
        neAvalVec.clear();
    }
    
    while (hIbf->GetMaximum() < 100) hIbf->Rebin(2);
    hIbf->Scale(1/hIbf->GetMaximum());
    hIbf->SetMaximum(1.35);
    TF1* fIbf = GetFitCurve(hIbf);
    
    bool gaussian = (fIbf->GetParameter(0)) > 0;
    if (!gaussian) fIbf = GetFitCurve(hIbf,false);
    
    hIbf->SetTitle("IBF");
    hIbf->SetLineColor(kBlue);
    hIbf->Draw("hist same");
    fIbf->SetLineColor(kBlue);
    fIbf->Draw("same");
    
    // Using charges
    
    TTree* tChargeReadout = (TTree*)fSignal->Get(Form("tInducedCharge_%d", readoutElectrode));
    TTree* tChargeDrift = (TTree*)fSignal->Get(Form("tInducedCharge_%d", driftElectrode));
    Double_t chargeDrift, chargeReadout, ionChargeDrift, ionChargeReadout;
    tChargeReadout->SetBranchAddress("totalInducedCharge", &chargeReadout);
    tChargeDrift->SetBranchAddress("totalInducedCharge", &chargeDrift);
    tChargeReadout->SetBranchAddress("ionInducedCharge", &ionChargeReadout);
    tChargeDrift->SetBranchAddress("ionInducedCharge", &ionChargeDrift);
    
    int nChargeReadout = tChargeReadout->GetEntries();
    int nChargeDrift = tChargeDrift->GetEntries();
    if (nChargeReadout != nChargeDrift) {
        std::cout << "nChargeReadout != nChargeDrift" << std::endl;
        return;
    }
    int nCharge = nChargeDrift;
    
    int nBins = 10000;
    TH1F* hIbfCharge = new TH1F("hIbfCharge", "hIbfCharge", nBins, 0, 100);
    TH1F* hIbfIonCharge = new TH1F("hIbfIonCharge", "hIbfIonCharge", nBins, 0, 100);
    for (int l = 0; l<nCharge; l++) {
        tChargeReadout->GetEntry(l);
        tChargeDrift->GetEntry(l);
        //tAvalanche->GetEntry(l);
        //if (nAmplification>1) {
        if (abs(chargeReadout)>10) {    // there is signal! it won't be 0 divided by 0
            // + if there's no signal, not relevant to compute an IBF
            hIbfCharge->Fill(abs(chargeDrift/chargeReadout*100.));
            hIbfIonCharge->Fill(abs(ionChargeDrift/ionChargeReadout*100.));
        }
    }
    while (hIbfCharge->GetMaximum() < 100) hIbfCharge->Rebin(2);
    while (hIbfIonCharge->GetMaximum() < 100) hIbfIonCharge->Rebin(2);
    hIbfCharge->Scale(1./hIbfCharge->GetMaximum());
    hIbfCharge->SetMaximum(1.35);
    hIbfIonCharge->Scale(1./hIbfIonCharge->GetMaximum());
    hIbfIonCharge->SetMaximum(1.35);
    TF1* fIbfCharge = GetFitCurve(hIbfCharge, gaussian);
    TF1* fIbfIonCharge = GetFitCurve(hIbfIonCharge, gaussian);
    
    // Change x axis
    /*
    hIbf->GetXaxis()->SetRangeUser(0, 10.);
    hIbfCharge->GetXaxis()->SetRangeUser(0, 10.);
    hIbfIonCharge->GetXaxis()->SetRangeUser(0, 10.);
    */
    hIbf->GetXaxis()->SetRangeUser(0, 2*fIbfCharge->GetParameter(0));
    hIbfCharge->GetXaxis()->SetRangeUser(0, 2*fIbfCharge->GetParameter(0));
    hIbfIonCharge->GetXaxis()->SetRangeUser(0, 2*fIbfIonCharge->GetParameter(0));
    if (fIbf->GetParameter(0) < 0.9) {
        hIbf->GetXaxis()->SetRangeUser(0, 3.);
        hIbfCharge->GetXaxis()->SetRangeUser(0, 3.);
        hIbfIonCharge->GetXaxis()->SetRangeUser(0, 3.);
    }
    if (fIbf->GetParameter(0) < 0.5) {
        hIbf->GetXaxis()->SetRangeUser(0, 1.);
        hIbfCharge->GetXaxis()->SetRangeUser(0, 1.);
        hIbfIonCharge->GetXaxis()->SetRangeUser(0, 1.);
    }
    if (fIbf->GetParameter(0) < 0.1) {
        hIbf->GetXaxis()->SetRangeUser(0, 0.2);
        hIbfCharge->GetXaxis()->SetRangeUser(0, 0.2);
        hIbfIonCharge->GetXaxis()->SetRangeUser(0, 0.2);
        hIbf->GetXaxis()->SetMaxDigits(3);
        hIbfCharge->GetXaxis()->SetMaxDigits(3);
        hIbfIonCharge->GetXaxis()->SetMaxDigits(3);
    }
    hIbf->SetXTitle("IBF (%)");
    hIbfCharge->SetXTitle("IBF (%)");
    hIbfIonCharge->SetXTitle("IBF (%)");
    
    hIbfCharge->SetLineColor(7);
    fIbfCharge->SetLineColor(7);
    
    hIbfIonCharge->SetLineColor(6);
    fIbfIonCharge->SetLineColor(6);
    
    hIbfCharge->Draw("hist same");
    fIbfCharge->Draw("same");
    
    hIbfIonCharge->Draw("hist same");
    fIbfIonCharge->Draw("same");
    
    // Finally legend
    
    TString txtIbf1 = Form("IBF = (%.2f #pm %.2f) %s", fIbf->GetParameter(0), fIbf->GetParError(0), "%");
    TString txtIbf2 = Form("IBF = (%.2f #pm %.2f) %s", fIbfCharge->GetParameter(0), fIbfCharge->GetParError(0), "%");
    TString txtIbf3 = Form("IBF = (%.2f #pm %.2f) %s", fIbfIonCharge->GetParameter(0), fIbfIonCharge->GetParError(0), "%");
    TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
    legend->SetMargin(0.15);
    legend->AddEntry(hIbf,"IBF ratio: " + txtIbf1,"l");
    legend->AddEntry(hIbfCharge,"All induced charges: " + txtIbf2,"l");
    legend->AddEntry(hIbfIonCharge,"Induced ion charges: " + txtIbf3,"l");
    legend->SetTextSize(0.04);
    legend->Draw("same");
    
    
}


void DrawConvolutedIbf(TString fConvolutedName="") {
    
    // First draw convoluted spectrum
    TFile* fConvoluted = TFile::Open(fConvolutedName, "READ");
    
    TH1F* hFeIbf = (TH1F*)fConvoluted->Get("hFeIbf");
    TH1F* hFeIbfTotalCharge = (TH1F*)fConvoluted->Get("hFeIbfTotalCharge");
    TH1F* hFeIbfIonCharge = (TH1F*)fConvoluted->Get("hFeIbfIonCharge");
    hFeIbf->Scale(1/hFeIbf->GetMaximum());
    hFeIbfTotalCharge->Scale(1/hFeIbfTotalCharge->GetMaximum());
    hFeIbfIonCharge->Scale(1/hFeIbfIonCharge->GetMaximum());
    
    TF1* fFeIbf = GetFitCurve(hFeIbf);
    TF1* fFeIbfTotalCharge = GetFitCurve(hFeIbfTotalCharge);
    TF1* fFeIbfIonCharge = GetFitCurve(hFeIbfIonCharge);
    /*
    hFeIbf->GetXaxis()->SetRangeUser(0, 10.);
    hFeIbfTotalCharge->GetXaxis()->SetRangeUser(0, 10.);
    hFeIbfIonCharge->GetXaxis()->SetRangeUser(0, 10.);
     */
    hFeIbf->GetXaxis()->SetRangeUser(0, 2*fFeIbf->GetParameter(0));
    hFeIbfTotalCharge->GetXaxis()->SetRangeUser(0, 2*fFeIbfTotalCharge->GetParameter(0));
    hFeIbfIonCharge->GetXaxis()->SetRangeUser(0, 2*fFeIbfIonCharge->GetParameter(0));
    if (fFeIbf->GetParameter(0) < 0.9) {
        hFeIbf->GetXaxis()->SetRangeUser(0, 3.);
        hFeIbfTotalCharge->GetXaxis()->SetRangeUser(0, 3.);
        hFeIbfIonCharge->GetXaxis()->SetRangeUser(0, 3.);
    }
    if (fFeIbf->GetParameter(0) < 0.5) {
        hFeIbf->GetXaxis()->SetRangeUser(0, 1.);
        hFeIbfTotalCharge->GetXaxis()->SetRangeUser(0, 1.);
        hFeIbfIonCharge->GetXaxis()->SetRangeUser(0, 1.);
    }
    if (fFeIbf->GetParameter(0) < 0.1) {
        hFeIbf->GetXaxis()->SetRangeUser(0, 0.2);
        hFeIbfTotalCharge->GetXaxis()->SetRangeUser(0, 0.2);
        hFeIbfIonCharge->GetXaxis()->SetRangeUser(0, 0.2);
        hFeIbf->GetXaxis()->SetMaxDigits(3);
        hFeIbfTotalCharge->GetXaxis()->SetMaxDigits(3);
        hFeIbfIonCharge->GetXaxis()->SetMaxDigits(3);
    }
    hFeIbf->SetXTitle("IBF (%)");
    hFeIbfTotalCharge->SetXTitle("IBF (%)");
    hFeIbfIonCharge->SetXTitle("IBF (%)");
    
    hFeIbf->SetLineColor(12);
    fFeIbf->SetLineColor(12);
    hFeIbfTotalCharge->SetLineColor(9);
    fFeIbfTotalCharge->SetLineColor(9);
    hFeIbfIonCharge->SetLineColor(8);
    fFeIbfIonCharge->SetLineColor(8);
    
    hFeIbf->SetMaximum(1.35);
    hFeIbfTotalCharge->SetMaximum(1.35);
    hFeIbfIonCharge->SetMaximum(1.35);
    
    // Draw convoluted IBF
    /*
     hFeIbf->Draw("hist same");
     fFeIbf->Draw("same");
     
     hFeIbfTotalCharge->Draw("hist same");
     fFeIbfTotalCharge->Draw("same");
     */
    
    hFeIbfIonCharge->Draw("hist same");
    fFeIbfIonCharge->Draw("same");
    
    TString txtIbf1 = Form("IBF = (%.2f #pm %.2f) %s", fFeIbf->GetParameter(0), fFeIbf->GetParError(0), "%");
    TString txtIbf2 = Form("IBF = (%.2f #pm %.2f) %s", fFeIbfTotalCharge->GetParameter(0), fFeIbfTotalCharge->GetParError(0), "%");
    TString txtIbf3 = Form("IBF = (%.2f #pm %.2f) %s", fFeIbfIonCharge->GetParameter(0), fFeIbfIonCharge->GetParError(0), "%");
    TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
    legend->SetMargin(0.15);
    /*
     legend->AddEntry(fFeIbf,"IBF ratio: " + txtIbf1,"l");
     legend->AddEntry(fFeIbfTotalCharge,"All induced charges: " + txtIbf2,"l");
     */
    legend->AddEntry(fFeIbfIonCharge,"Induced ion charges: " + txtIbf3,"l");
    legend->SetTextSize(0.04);
    legend->Draw("same");
    
}



int Ibf(int modelNum, std::string gasName, std::vector<int> hvList, bool useFeSource) {
    
    time_t t0 = time(NULL);
    gStyle->SetTitleFontSize(.06);
    gStyle->SetTitleSize(.06);
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetLabelSize(.04, "XY");
    gStyle->SetMarkerSize(0.3);
    gStyle->SetTextSize(0.05);
    
    int electrodeNum = GetElectrodeNum(modelNum);
    if ((int)hvList.size() != electrodeNum-1) {cout << "Wrong hv input" << endl; return 0;}
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    // input and output files
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    TString fConvolutedName = path + "fe-spectrum-convoluted";
    TString fSignalName = path;
    if (useFeSource) fSignalName += "fe-signal";
    else fSignalName += "signal";
    //TString fSignalName = path + "signal";
    TString outputName = Form("Figures/model%d/ibf-%s", modelNum, gasName.c_str());
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
    DrawIbf(modelNum, fSignalName);
    
    cv->cd(2);
    DrawConvolutedIbf(fConvolutedName);
    
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


