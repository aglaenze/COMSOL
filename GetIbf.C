#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "Include/Utils.C"
#include "Include/Data.C"
#include "Include/Functions.C"



int GetIbf() {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    bool convoluted = true;     // to draw convoluted ibf with Fe source instead of "pure" ibf
    //____________________
    
    time_t t0 = time(NULL);
    
    gStyle->SetTitleFontSize(.06);
    gStyle->SetTitleSize(.06);
    gStyle->SetTextSize(.05);
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.3);
    gStyle->SetMarkerStyle(20);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetTextSize(.05);
    
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    
    // Get number of files to look at
    //Int_t num = 4;
    Int_t num = GetNumberOfFiles(path, "signal");
    Int_t num2 = int(num/2.);
    if (num/2.> num2) num2+=1;
    TCanvas* c2 = new TCanvas("c2", "c2", 400, 300*num2);
    c2->Divide(2, num2);
    
    Double_t ibfList[num], ibfErrorList[num], ibfList2[num], ibfErrorList2[num], ibfList3[num];
    Double_t ibfConvolutedList[num], ibfConvolutedErrorList[num], ibfConvolutedList2[num], ibfConvolutedErrorList2[num], ibfConvolutedList3[num], ibfConvolutedErrorList3[num];
    Double_t hvMeshList[num], hvDriftList[num], hvRatioList[num];
    Double_t fieldList[num];
    double damp = 0.0128, ddrift = 0.5;
    
    const Int_t dataNum = dataQuantity(gasName);
    
    int readoutElectrode = 0;
    int driftElectrode = 0;
    std::map <std::string, int, NoSorting> electrode;
    LoadElectrodeMap(modelNum, electrode);
    
    map<string, int>::iterator it = electrode.begin();
    for (it=electrode.begin(); it!=electrode.end(); ++it) {
        //cout << it->first << " => " << it->second << '\n';
        if (it->first == "pad") readoutElectrode = it->second;  // could be mesh, it depends on where you want to read
        else if (it->first == "drift") driftElectrode = it->second;
    }
    if (readoutElectrode == 0 || driftElectrode == 0) {cout << "Did not find drift or pad electrode" << endl; return 0;}
    
    for (int k = 0; k<num; k++) {
        int hvMesh = 340+20*k;
        int hvDrift = 200+hvMesh;
        //hvRatioList[k] = (double)hvMesh/(hvDrift-hvMesh)*(ddrift-damp)/damp;
        hvRatioList[k] = (double)hvMesh/(hvDrift-hvMesh)*ddrift/damp;
        
        TFile* fSignal = TFile::Open(Form("rootFiles/Ar-iC4H10/model%d/signal-%d-%d.root", modelNum, hvMesh, hvDrift), "READ");
        TFile* fConvoluted = TFile::Open(Form("rootFiles/Ar-iC4H10/model%d/fe-spectrum-convoluted-%d-%d.root", modelNum, hvMesh, hvDrift), "READ");
        
        // 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
        // 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
        // 3e étape: faire une ligne verticale de l'estimation du ratio des intégrales
        // 4e étape : dessiner la courbe ibf = f(field ratio)
        
        // 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
        TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
        int nAvalanche = tAvalanche->GetEntries();
        //Int_t nAmplification;
        vector <Int_t> *nAmpPtr = nullptr;
        vector <Int_t> nAmp = {};
        Int_t ni = 0, ionBackNum = 0;
        //tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
        tAvalanche->SetBranchAddress("amplificationElectrons", &nAmpPtr);
        tAvalanche->SetBranchAddress("ionNum", &ni);
        tAvalanche->SetBranchAddress("ionBackNum", &ionBackNum);
        
        TH1F* hIbf = new TH1F("ibf", "ibf", 2000, 0, 100);
        hIbf->SetXTitle("IBF (%)");
        for (int l = 0; l<nAvalanche; l++) {
            tAvalanche->GetEntry(l);
            nAmp = *nAmpPtr;
            if (nAmp[0]>1) hIbf->Fill((double)ionBackNum/ni*100.);
        }
        
        hIbf->Scale(1/hIbf->GetMaximum());
        hIbf->SetMaximum(1.2);
        TF1* f = GetFitIbf(hIbf);
        
        
        Int_t iBinMax = hIbf->GetMaximumBin();
        Double_t xMax = hIbf->GetXaxis()->GetBinCenter( iBinMax );
        
        double xPlotMax = 3*f->GetParameter(0);
        hIbf->GetXaxis()->SetRangeUser(0, xPlotMax);
        hvMeshList[k] = hvMesh;
        fieldList[k] = hvMesh/damp/1000.;
        ibfList[k] = f->GetParameter(0);
        ibfErrorList[k] = f->GetParError(0);
        
        // 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
        
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
            return 0;
        }
        int nCharge = nChargeDrift;
        
        TH1F* hIbfCharge = new TH1F("hIbfCharge", "hIbfCharge", 1000, 0, 100);
        TH1F* hIbfIonCharge = new TH1F("hIbfIonCharge", "hIbfIonCharge", 1000, 0, 100);
        for (int l = 0; l<nCharge; l++) {
            tChargeReadout->GetEntry(l);
            tChargeDrift->GetEntry(l);
            hIbfCharge->Fill(abs(chargeDrift/chargeReadout*100.));
            hIbfIonCharge->Fill(abs(ionChargeDrift/ionChargeReadout*100.));
        }
        hIbfCharge->Scale(1/hIbfCharge->GetMaximum());
        hIbfCharge->SetMaximum(1.2);
        hIbfIonCharge->Scale(1/hIbfIonCharge->GetMaximum());
        TF1* fIbfCharge = GetFitIbf(hIbfCharge);
        TF1* fIbfIonCharge = GetFitIbf(hIbfIonCharge, true);
        
        
        /*
         Int_t iBinMax2 = hIbfCharge->GetMaximumBin();
         Double_t xMax2 = hIbfCharge->GetXaxis()->GetBinCenter( iBinMax );
         hIbfCharge->GetXaxis()->SetRangeUser(0, xMax + 3*hIbf->GetRMS());
         */
        
        hIbfCharge->GetXaxis()->SetRangeUser(0, xPlotMax);
        ibfList2[k] = fIbfCharge->GetParameter(0);
        ibfErrorList2[k] = fIbfCharge->GetParError(0);
        
        // 3e étape (nouvelle étape intermédiaire) : récupérer les histos IBF convolués, et les fitter
        TH1F* hFeIbf = (TH1F*)fConvoluted->Get("hFeIbf");
        TH1F* hFeIbfTotalCharge = (TH1F*)fConvoluted->Get("hFeIbfTotalCharge");
        TH1F* hFeIbfIonCharge = (TH1F*)fConvoluted->Get("hFeIbfIonCharge");
        hFeIbf->Scale(1/hFeIbf->GetMaximum());
        hFeIbfTotalCharge->Scale(1/hFeIbfTotalCharge->GetMaximum());
        hFeIbfIonCharge->Scale(1/hFeIbfIonCharge->GetMaximum());
        hFeIbf->SetMaximum(1.2);
        hFeIbfTotalCharge->SetMaximum(1.2);
        hFeIbfIonCharge->SetMaximum(1.2);
        
        if (convoluted) {
            iBinMax = hFeIbf->GetMaximumBin();
            xMax = hFeIbf->GetXaxis()->GetBinCenter( iBinMax );
        }
        
        TF1* fFeIbf = GetFitIbf(hFeIbf);
        TF1* fFeIbfTotalCharge = GetFitIbf(hFeIbfTotalCharge);
        TF1* fFeIbfIonCharge = GetFitIbf(hFeIbfIonCharge);
        if (convoluted) xPlotMax = 3*fFeIbf->GetParameter(0);
        
        hFeIbf->GetXaxis()->SetRangeUser(0, xPlotMax);
        hFeIbfIonCharge->GetXaxis()->SetRangeUser(0, xPlotMax);
        hFeIbfTotalCharge->GetXaxis()->SetRangeUser(0, xPlotMax);
        
        // Extract IBF values
        ibfConvolutedList[k] = fFeIbf->GetParameter(0);
        ibfConvolutedErrorList[k] = fFeIbf->GetParError(0);
        
        ibfConvolutedList2[k] = fFeIbfTotalCharge->GetParameter(0);
        ibfConvolutedErrorList2[k] = fFeIbfTotalCharge->GetParError(0);
        
        ibfConvolutedList3[k] = fFeIbfIonCharge->GetParameter(0);
        ibfConvolutedErrorList3[k] = fFeIbfIonCharge->GetParError(0);
        
        // 4e étape : dessiner les courbes ibf = f(field ratio)
        c2->cd(k+1);
        // Upper plot will be in pad1
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 0.95);
        pad1->SetTopMargin(0.2);
        pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->SetLeftMargin(0.15);
        //pad1->SetGridx();         // Vertical grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        if (convoluted) {
            hFeIbf->Draw("hist");
            fFeIbf->Draw("same");
        }
        else {
            hIbf->Draw("hist");
            f->Draw("same");
        }
        
        
        // Add text to frame
        double xPos = 1.6*xMax;
        TLatex* txt = new TLatex(xPos,0.6,Form("#bf{IBF = %.2f #pm %.2f %s}", f->GetParameter(0), f->GetParError(0), "%"));
        TLatex* txtConv = new TLatex(xPos,0.6,Form("#bf{IBF = %.2f #pm %.2f %s}",fFeIbf->GetParameter(0), fFeIbf->GetParError(0), "%"));
        TLatex* txt2 = new TLatex(xPos,0.5,Form("#bf{Field ratio = %.2f}", hvRatioList[k]));
        txt->SetTextSize(0.05) ;
        txt2->SetTextSize(0.05) ;
        if (convoluted) txtConv->Draw();
        else txt->Draw();
        txt2->Draw();
        
        // lower plot will be in pad
        c2->cd(k+1);          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.5);
        pad2->SetTopMargin(0);
        pad2->SetLeftMargin(0.15);
        pad2->SetBottomMargin(0.2);
        //pad2->SetGridx(); // vertical grid
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad
        if (convoluted) {
            hFeIbfTotalCharge->Draw("hist");
            hFeIbfIonCharge->SetLineColor(6);
            hFeIbfIonCharge->Draw("hist same");
            fFeIbfTotalCharge->Draw("same");
        }
        else {
            hIbfCharge->Draw("hist");
            hIbfIonCharge->SetLineColor(6);
            hIbfIonCharge->Draw("hist same");
            fIbfCharge->Draw("same");
        }
        
        TLegend* legend = new TLegend(0.6,0.7,0.9,0.9);
        legend->SetTextSize(0.045);
        if (convoluted) {
            legend->AddEntry(hFeIbfTotalCharge,"All induced charges","l");
            legend->AddEntry(hFeIbfIonCharge,"Induced ion charges","l");
        }
        else {
            legend->AddEntry(hIbfCharge,"All induced charges","l");
            legend->AddEntry(hIbfIonCharge,"Induced ion charges","l");
        }
        legend->Draw("same");
        
        TLatex* txt3 = new TLatex(xPos, 0.6, Form("#bf{IBF = %.2f #pm %.2f %s}", fIbfCharge->GetParameter(0), fIbfCharge->GetParError(0), "%"));
        TLatex* txt3Conv = new TLatex(xPos, 0.6, Form("#bf{IBF = %.2f #pm %.2f %s}", fFeIbfTotalCharge->GetParameter(0), fFeIbfTotalCharge->GetParError(0), "%"));
        TLatex* txt4 = new TLatex(xPos, 0.5, Form("#bf{Field ratio = %.2f}", hvRatioList[k]));
        txt3->SetTextSize(0.05) ;
        txt4->SetTextSize(0.05) ;
        if (convoluted) txt3Conv->Draw();
        else txt3->Draw();
        txt4->Draw();
        
        
        c2->cd(k+1);
        //ibfLine->Draw("same");
    }
    c2->SaveAs(Form("Figures/model%d/ibf_%s.pdf", modelNum, gasName.c_str()));
    
    
    TCanvas* c5 = new TCanvas("c5", "c5", 600, 300);
    c5->Divide(2);
    c5->SetGrid();
    gPad->SetGrid();
    
    c5->cd(1);
    gPad->SetGrid();
    
    TGraphErrors* grSim1 = new TGraphErrors(num, fieldList, ibfList, 0, ibfErrorList);
    grSim1->SetTitle("IBF curve in the Micromegas");
    grSim1->GetXaxis()->SetTitle( "E_{amp} (kV/cm)" );
    grSim1->GetYaxis()->SetTitle( "IBF (%)" );
    //grSim1->GetXaxis()->SetLimits(0, 4);  // along X axis
    grSim1->GetHistogram()->SetMinimum(0.);   // along Y axis
    grSim1->GetHistogram()->SetMaximum(5.);   // along Y axis
    grSim1->SetMarkerColor(1);
    //grSim1->GetXaxis()->SetLimits(hvMeshList[0]-5, hvMeshList[num]+5);
    //grSim1->Draw("ALP");
    
    TGraphErrors* grSim2 = new TGraphErrors(num, fieldList, ibfList2, 0, ibfErrorList2);
    grSim2->SetMarkerColor(2);
    //grSim2->Draw("LP same");
    
    TGraphErrors* grSimConvoluted = new TGraphErrors(num, fieldList, ibfConvolutedList, 0, ibfConvolutedErrorList);
    //grSimConvoluted->SetMarkerColor(4);
    grSimConvoluted->SetTitle("IBF curve in the Micromegas");
    grSimConvoluted->GetXaxis()->SetTitle( "E_{amp} (kV/cm)" );
    grSimConvoluted->GetYaxis()->SetTitle( "IBF (%)" );
    grSimConvoluted->GetHistogram()->SetMinimum(0.);   // along Y axis
    grSimConvoluted->GetHistogram()->SetMaximum(3.5);   // along Y axis
    grSimConvoluted->Draw("AP same");
    
    //TGraphErrors* grSimConvoluted2 = new TGraphErrors(num, fieldList, ibfConvolutedList2, 0, ibfConvolutedErrorList2);
    TGraphErrors* grSimConvoluted2 = new TGraphErrors(num, fieldList, ibfList, 0, ibfErrorList);
    //grSimConvoluted2->SetMarkerColor(7);
    //grSimConvoluted2->Draw("LP same");
    grSimConvoluted2->SetTitle("IBF curve in the Micromegas");
    grSimConvoluted2->GetXaxis()->SetTitle( "E_{amp} (kV/cm)" );
    grSimConvoluted2->GetYaxis()->SetTitle( "IBF (%)" );
    grSimConvoluted2->GetHistogram()->SetMinimum(0.);   // along Y axis
    grSimConvoluted2->GetHistogram()->SetMaximum(3.5);   // along Y axis
    //grSimConvoluted2->Draw("AP same");
    
    TGraphErrors* grSimConvoluted3 = new TGraphErrors(num, fieldList, ibfConvolutedList3, 0, ibfConvolutedErrorList3);
    grSimConvoluted3->SetMarkerColor(6);
    //grSimConvoluted3->Draw("LP same");
    
    // Same with data
    //const Int_t dataNum = dataQuantity(gasName);
    Double_t hvMeshListIBF1[dataNum], hvDriftListIBF1[dataNum], ionBackFlowCorrectedVect1[dataNum], ionBackFlowCorrectedErrorVect1[dataNum], ionBackFlowCorrectedVect1_old[dataNum], ionBackFlowCorrectedErrorVect1_old[dataNum];
    LoadIbfData(gasName, dataNum, hvMeshListIBF1, hvDriftListIBF1, ionBackFlowCorrectedVect1, ionBackFlowCorrectedErrorVect1, ionBackFlowCorrectedVect1_old, ionBackFlowCorrectedErrorVect1_old);
    
    Double_t ratioListData[dataNum], fieldListData[dataNum];
    for (int i = 0; i < dataNum; i++) {
        ratioListData[i] = hvMeshListIBF1[i]/(hvDriftListIBF1[i]-hvMeshListIBF1[i]) * (ddrift-damp)/damp;
        fieldListData[i] = hvMeshListIBF1[i]/damp;
    }
    
    TGraphErrors* grData1 = new TGraphErrors(dataNum, fieldListData, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
    grData1->SetMarkerColor(4);
    //grData1->Draw("LP same");
    TGraphErrors* grData2 = new TGraphErrors(dataNum, fieldListData, ionBackFlowCorrectedVect1_old, 0, ionBackFlowCorrectedErrorVect1_old );
    grData2->SetMarkerColor(5);
    //grData2->Draw("LP same");
    
    TLegend* legend = new TLegend(0.45,0.65,0.9,0.9);
    legend->SetTextSize(0.03);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    /*legend->AddEntry(grData1,"Data","lp");
     legend->AddEntry(grData2,"Data old","lp");
     */
    legend->AddEntry(grData2,"Data","lp");
    legend->AddEntry(grSim1,"Simulation (counting IBF)","lp");
    legend->AddEntry(grSim2, "Simulation (induced charge)", "lp");
    legend->AddEntry(grSimConvoluted, "Convolution (counting IBF)", "lp");
    legend->AddEntry(grSimConvoluted2, "Convolution (total induced charge)", "lp");
    legend->AddEntry(grSimConvoluted3, "Convolution (induced ion charge)", "lp");
    //legend->Draw();
    
    
    c5->cd(2);
    gPad->SetGrid();
    TGraphErrors* grSimRatio1 = nullptr;
    //TGraphErrors* grSimRatio1 = new TGraphErrors(num, hvRatioList, ibfConvolutedList2, 0, ibfConvolutedErrorList2);
    if (convoluted) grSimRatio1 = new TGraphErrors(num, hvRatioList, ibfConvolutedList, 0, ibfConvolutedErrorList);
    else grSimRatio1 = new TGraphErrors(num, hvRatioList, ibfList, 0, ibfErrorList);
    grSimRatio1->SetTitle("IBF = f(E_{amp}/E_{drift})");
    grSimRatio1->GetXaxis()->SetTitle( "E_{amp}/E_{drift}" );
    grSimRatio1->GetYaxis()->SetTitle( "IBF (%)" );
    //grSimRatio1->GetXaxis()->SetLimits(0, 4);  // along X axis
    grSimRatio1->GetHistogram()->SetMinimum(0.);   // along Y axis
    grSimRatio1->GetHistogram()->SetMaximum(3.5);   // along Y axis
    //grSimRatio1->SetMarkerColor(1);
    //grSimRatio1->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num]+5);
    grSimRatio1->Draw("AP");
    
    TGraphErrors* grSimRatio2 = new TGraphErrors(num, hvRatioList, ibfList2, 0, ibfErrorList);
    grSimRatio2->SetTitle("IBF = f(E ratio)");
    grSimRatio2->GetXaxis()->SetTitle( "E_{amp}/E_{drift}" );
    grSimRatio2->GetYaxis()->SetTitle( "IBF (%)" );
    //grSimRatio2->GetXaxis()->SetLimits(0, 4);  // along X axis
    grSimRatio2->GetHistogram()->SetMinimum(0.);   // along Y axis
    grSimRatio2->GetHistogram()->SetMaximum(4.);   // along Y axis
    grSimRatio2->SetMarkerColor(2);
    //grSimRatio2->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num]+5);
    //grSimRatio2->Draw("LP same");
    //grSimRatio2->Draw("AP");
    
    
    TGraphErrors* grDataRatio1 = new TGraphErrors(dataNum, ratioListData, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
    grDataRatio1->SetMarkerStyle(20);
    grDataRatio1->SetMarkerColor(4);
    //grDataRatio1->Draw("LP same");
    TGraphErrors* grDataRatio2 = new TGraphErrors(dataNum, ratioListData, ionBackFlowCorrectedVect1_old, 0, ionBackFlowCorrectedErrorVect1_old );
    grDataRatio2->SetMarkerStyle(20);
    grDataRatio2->SetMarkerColor(5);
    //grDataRatio2->Draw("LP same");
    
    TLegend* legendRatio = new TLegend(0.45,0.65,0.9,0.9);
    legendRatio->SetTextSize(0.03);
    //legendRatio->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    /*legendRatio->AddEntry(grDataRatio1, "Data", "lp");
     legendRatio->AddEntry(grDataRatio2, "Data old", "lp");
     */
    legendRatio->AddEntry(grDataRatio2, "Data", "lp");
    legendRatio->AddEntry(grSimRatio1, "Simulation (counting IBF)", "lp");
    legendRatio->AddEntry(grSimRatio2, "Simulation (induced charge)", "lp");
    //legendRatio->Draw();
    
    c5->SaveAs(Form("Figures/model%d/IBFCurve-%s.pdf", modelNum, gasName.c_str()));
    
    
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    vector<double> ibfVec = {};
    vector<double> ibfErrorVec = {};
    vector<double> fieldVec = {};
    vector<double> fieldRatioVec = {};
    for (int k = 0; k < num; k++) {
        if (convoluted) {
            ibfVec.push_back(ibfConvolutedList[k]);
            ibfErrorVec.push_back(ibfConvolutedErrorList[k]);
        }
        else {
            ibfVec.push_back(ibfList[k]);
            ibfErrorVec.push_back(ibfErrorList[k]);
        }
        fieldVec.push_back(fieldList[k]);
        fieldRatioVec.push_back(hvRatioList[k]);
    }
    PrintListVec( "ibfListSim", ibfVec, "%.3f" );
    PrintListVec( "ibfErrorListSim", ibfErrorVec, "%.3f" );
    PrintListVec( "fieldListSim", fieldVec, "%.2f" );
    PrintListVec( "fieldRatioListSim", fieldRatioVec, "%.2f" );
    
    
    return 0;
}


