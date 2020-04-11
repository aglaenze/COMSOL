#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "_Utils.C"
#include "_Data.C"



Double_t FitGauss( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
    return  par[2]*TMath::Gaus( x[0], par[0], par[1]); }

TF1* GetFitIbf(TH1F* h) {
    Int_t iBinMax = h->GetMaximumBin();
    Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
    
    std::cout << "xMax = " << xMax << std::endl;
    std::cout << "maximum = " << h->GetMaximum() << std::endl;
    
    Int_t fitRangeMin = xMax - 0.6 * h->GetRMS();
    Int_t fitRangeMax = xMax + 0.6 * h->GetRMS();

    TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
    f->SetParNames("Mean", "Sigma", "Amplitude");
    f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
    
    h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
    return f;
}


int GetIbf() {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    //____________________

    time_t t0 = time(NULL);
    
    gStyle->SetTitleFontSize(.06);
    gStyle->SetTitleSize(.06);

     gStyle->SetOptStat(0);
     gStyle->SetTitleFontSize(.05);
     gStyle->SetTitleXSize(.05);
     gStyle->SetTitleYSize(.05);
     gStyle->SetLabelSize(.05, "XY");

    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);

    // Get number of files to look at
    Int_t num = 1;
    /*
    Int_t num = GetNumberOfFiles(path, "signal");
    Int_t num2 = int(num/2.);
    if (num/2.> num2) num2+=1;
     */
    
    Double_t ibfList[num], ibfErrorList[num], ibfList2[num], ibfErrorList2[num];
    Double_t hvMeshList[num], hvDriftList[num], hvRatioList[num];
    double damp = 0.0128, ddrift = 0.5;
    
    int k =0;
    int hvMesh = 340;
    int hvDrift = 200+hvMesh;
    
    TFile* fSignal = TFile::Open(Form("rootFiles/Ar-iC4H10/model1/signal-%d-%d.root", hvMesh, hvDrift), "READ");
    
    // 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
    // 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
    // 3e étape: faire une ligne verticale de l'estimation du ratio des intégrales
    // 4e étape : dessiner la courbe ibf = f(field ratio)
    
    // 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    Double_t ibfRatio;
    tAvalanche->SetBranchAddress("ibfRatio", &ibfRatio);
    int nAvalanche = tAvalanche->GetEntries();
    
    TH1F* hIbf = new TH1F("ibf", "ibf", 1000, 0, 100);
    hIbf->SetXTitle("IBF (%)");
    for (int l = 0; l<nAvalanche; l++) {
        tAvalanche->GetEntry(l);
        hIbf->Fill(ibfRatio*100.);
    }
    
    TF1* f = GetFitIbf(hIbf);
       
    /*
    Int_t iBinMax = hIbf->GetMaximumBin();
    Double_t xMax = hIbf->GetXaxis()->GetBinCenter( iBinMax );
    hIbf->GetXaxis()->SetRangeUser(0, xMax + 3*hIbf->GetRMS());
     */
    hIbf->GetXaxis()->SetRangeUser(0, 10.);
    hvMeshList[k] = hvMesh;
    ibfList[k] = f->GetParameter(0);
    ibfErrorList[k] = f->GetParError(0);

    
    // 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
    std::map <std::string, int> electrode;
    LoadElectrodeMap(modelNum, electrode);
    
    int readoutElectrode = electrode["pad"];    // could be mesh, it depends on wher you want to read
    int driftElectrode = electrode["drift"];
    
    TTree* tChargeReadout = (TTree*)fSignal->Get(Form("tInducedCharge_%d", readoutElectrode));
    TTree* tChargeDrift = (TTree*)fSignal->Get(Form("tInducedCharge_%d", driftElectrode));
    Double_t chargeDrift, chargeReadout;
    tChargeReadout->SetBranchAddress("totalInducedCharge", &chargeReadout);
    tChargeDrift->SetBranchAddress("totalInducedCharge", &chargeDrift);
    
    int nChargeReadout = tChargeReadout->GetEntries();
    int nChargeDrift = tChargeDrift->GetEntries();
    if (nChargeReadout != nChargeDrift) {
        std::cout << "nChargeReadout != nChargeDrift" << std::endl;
        return 0;
    }
    int nCharge = nChargeDrift;
    
    TH1F* hIbfCharge = new TH1F("hIbfCharge", "hIbfCharge", 1000, 0, 100);
    for (int l = 0; l<nCharge; l++) {
        tChargeReadout->GetEntry(l);
        tChargeDrift->GetEntry(l);
        hIbfCharge->Fill(abs(chargeDrift/chargeReadout*100.));
    }
    TF1* f2 = GetFitIbf(hIbfCharge);
    

    /*
    Int_t iBinMax2 = hIbfCharge->GetMaximumBin();
    Double_t xMax = hIbfCharge->GetXaxis()->GetBinCenter( iBinMax );
    hIbfCharge->GetXaxis()->SetRangeUser(0, xMax + 3*hIbf->GetRMS());
     */
    hIbfCharge->GetXaxis()->SetRangeUser(0, 10.);
    ibfList2[k] = f2->GetParameter(0);
    ibfErrorList2[k] = f2->GetParError(0);
    
    // 3e étape: dessiner une ligne verticale de l'estimation du ratio des intégrales
    
    TTree* tSignalReadout = (TTree*)fSignal->Get(Form("tSignal_%d", readoutElectrode));
    TTree* tSignalDrift = (TTree*)fSignal->Get(Form("tSignal_%d", driftElectrode));
    Double_t currentDriftIntMax = tSignalDrift->GetMaximum("totalCurrentInt");
    Double_t currentReadoutIntMax = tSignalReadout->GetMaximum("totalCurrentInt");
    if (currentReadoutIntMax<0) currentReadoutIntMax = abs(tSignalReadout->GetMinimum("totalCurrentInt"));
                                               
    Double_t ibfEstimation = 100*currentDriftIntMax/currentReadoutIntMax;
    
    Double_t yMax = 10;
    TLine* ibfLine = new TLine(ibfEstimation, 0, ibfEstimation, yMax);
    ibfLine->SetLineColor(2);
    

    // 4e étape : dessiner les courbes ibf = f(field ratio)
    TCanvas* c2 = new TCanvas("c2", "c2", 200, 300);
    // Now draw both spectra
    c2->cd(k+1);
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 0.95);
    pad1->SetTopMargin(0.2);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetLeftMargin(0.15);
    //pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    hIbf->Draw("hist");
    f->Draw("same");
    /*
    // Add text to frame
    TLatex* txt = new TLatex(0.5*xMax,1,Form("Gain = %.3g #pm %.3f", f->GetParameter(0), f->GetParError(0)));
    TLatex* txt2 = new TLatex(0.5*xMax,0.9,Form("Field ratio = %.2f", (double)hvMesh/(hvDrift-hvMesh)));
    txt->SetTextSize(0.05) ;
    txt2->SetTextSize(0.05) ;
    txt->Draw();
    txt2->Draw();
     */
        
    // lower plot will be in pad
    c2->cd(k+1);          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.5);
    pad2->SetTopMargin(0);
    pad2->SetLeftMargin(0.15);
    pad2->SetBottomMargin(0.2);
    //pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    hIbfCharge->Draw("hist");
    f2->Draw("same");
    /*
    TLatex* txt3 = new TLatex(0.5*xMax,1,Form("Gain = %.3g #pm %.3f", f2->GetParameter(0), f2->GetParError(0)));
    TLatex* txt4 = new TLatex(0.5*xMax,0.9,Form("Field ratio = %.2f", (double)hvMesh/(hvDrift-hvMesh)));
    txt3->SetTextSize(0.05) ;
    txt4->SetTextSize(0.05) ;
    txt3->Draw();
    txt4->Draw();
     */

    c2->cd(k+1);
    ibfLine->Draw("same");
    c2->SaveAs(Form("Figures/ibf_%s_model%d.pdf", gasName.c_str(), modelNum));
    return 0;

    
        TCanvas* c5 = new TCanvas("c5", "c5", 200, 500);
        c5->SetGrid();
        c5->Divide(1,3);
        
        c5->cd(1);
        TGraphErrors* gr = new TGraphErrors(num, hvMeshList, ibfList, 0, ibfErrorList);
        gr->SetTitle("IBF curve in the Micromegas");
        gr->GetXaxis()->SetTitle( "V_{mesh}" );
        gr->GetYaxis()->SetTitle( "IBF (%)" );
        //gr->GetXaxis()->SetLimits(0, 4);  // along X axis
        gr->GetHistogram()->SetMinimum(0.);   // along Y axis
        gr->GetHistogram()->SetMaximum(4.);   // along Y axis
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.3);
        //gr->GetXaxis()->SetLimits(hvMeshList[0]-5, hvMeshList[num]+5);
        gr->Draw("ACP");
        
        TGraph* tgIBFsignal = new TGraph(numS, hvMeshListS, ibfSignalList );
        tgIBFsignal->SetMarkerStyle(20);
        tgIBFsignal->SetMarkerSize(0.3);
        tgIBFsignal->SetMarkerColor(3);
        //tgIBFsignal->Draw("CP same");
        
        // Same with data
        //const Int_t dataNum = dataQuantity(gasName);
        Double_t hvMeshListIBF1[dataNum], hvDriftListIBF1[dataNum], ionBackFlowCorrectedVect1[dataNum], ionBackFlowCorrectedErrorVect1[dataNum], ionBackFlowCorrectedVect1_old[dataNum], ionBackFlowCorrectedErrorVect1_old[dataNum];
        LoadIbfData(gasName, dataNum, hvMeshListIBF1, hvDriftListIBF1, ionBackFlowCorrectedVect1, ionBackFlowCorrectedErrorVect1, ionBackFlowCorrectedVect1_old, ionBackFlowCorrectedErrorVect1_old);
        //Double_t hvMeshListData[dataNum], hvDriftListData[dataNum], gainListData[dataNum], gainErrorListData_old[dataNum], gainListData_old[dataNum], gainErrorListData[dataNum];
        //LoadGainData(gasName, dataNum, hvMeshListData, hvDriftListData, gainListData, gainErrorListData, gainListData_old, gainErrorListData_old);
        
        
        Double_t ratioListData[dataNum];
        for (int i = 0; i < dataNum; i++) {ratioListData[i] = hvMeshListIBF1[i]/(hvDriftListIBF1[i]-hvMeshListIBF1[i]) * (ddrift-damp)/damp;}
         
        TGraphErrors* tgIBF = new TGraphErrors(dataNum, hvMeshListIBF1, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
        tgIBF->SetMarkerStyle(20);
        tgIBF->SetMarkerSize(0.3);
        tgIBF->SetMarkerColor(2);
        tgIBF->Draw("CP same");
        TGraphErrors* tgIBF2 = new TGraphErrors(dataNum, hvMeshListIBF1, ionBackFlowCorrectedVect1_old, 0, ionBackFlowCorrectedErrorVect1_old );
        tgIBF2->SetMarkerStyle(20);
        tgIBF2->SetMarkerSize(0.3);
        tgIBF2->SetMarkerColor(3);
        tgIBF2->Draw("CP same");
     
        TLegend* legend = new TLegend(0.5,0.65,0.9,0.9);
        //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend->AddEntry(tgIBF,"Data","lp");
        legend->AddEntry(tgIBF2,"Data old","lp");
        legend->AddEntry(gr,"Simulation (counting IBF)","lp");
        //legend->AddEntry(tgIBFsignal, "Simulation (current signals)", "lp");
        legend->Draw();
       

        c5->cd(2);
        TGraphErrors* gr2 = new TGraphErrors(num, hvRatioList, ibfList, 0, ibfErrorList);
           gr2->SetTitle("IBF = f(E ratio)");
           gr2->GetXaxis()->SetTitle( "E_{amp}/E_{drift}" );
           gr2->GetYaxis()->SetTitle( "IBF (%)" );
           //gr2->GetXaxis()->SetLimits(0, 4);  // along X axis
           gr2->GetHistogram()->SetMinimum(0.);   // along Y axis
           gr2->GetHistogram()->SetMaximum(4.);   // along Y axis
           gr2->SetMarkerStyle(20);
        gr2->SetMarkerSize(0.3);
           //gr2->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num]+5);
           gr2->Draw("ACP");


           TGraphErrors* tgIBFr = new TGraphErrors(dataNum, ratioListData, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
           tgIBFr->SetMarkerStyle(20);
        tgIBFr->SetMarkerSize(0.3);
           tgIBFr->SetMarkerColor(2);
           tgIBFr->Draw("CP same");
           TGraphErrors* tgIBF2r = new TGraphErrors(dataNum, ratioListData, ionBackFlowCorrectedVect1_old, 0, ionBackFlowCorrectedErrorVect1_old );
           tgIBF2r->SetMarkerStyle(20);
        tgIBF2r->SetMarkerSize(0.3);
           tgIBF2r->SetMarkerColor(3);
           tgIBF2r->Draw("CP same");
        
           TLegend* legendr = new TLegend(0.5,0.65,0.9,0.9);
           //legendr->SetHeader("The Legend Title","C"); // option "C" allows to center the header
           legendr->AddEntry(tgIBF, "Data", "lp");
           legendr->AddEntry(tgIBF2, "Data old", "lp");
           legendr->AddEntry(gr,"Simulation (counting IBF)", "lp");
           //legendr->AddEntry(tgIBFsignal, "Simulation (current signals)", "lp");
           legendr->Draw();
        
        c5->cd(3);
        gPad->SetLogx();
        
        TGraphErrors* gr3 = new TGraphErrors(num, gainList, ibfList, gainErrorList, ibfErrorList);
           gr3->SetTitle("IBF = f(Gain)");
           gr3->GetXaxis()->SetTitle( "Gain" );
           gr3->GetYaxis()->SetTitle( "IBF (%)" );
           //gr3->GetXaxis()->SetLimits(0, 4);  // along X axis
           gr3->GetHistogram()->SetMinimum(0.);   // along Y axis
           gr3->GetHistogram()->SetMaximum(4.);   // along Y axis
           gr3->SetMarkerStyle(20);
        gr3->SetMarkerSize(0.3);
           //gr3->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num]+5);
           gr3->Draw("ACP");


           TGraphErrors* tgIBFg = new TGraphErrors(dataNum, gainListData, ionBackFlowCorrectedVect1, gainErrorListData, ionBackFlowCorrectedErrorVect1 );
           tgIBFg->SetMarkerStyle(20);
        tgIBFg->SetMarkerSize(0.3);
           tgIBFg->SetMarkerColor(2);
           tgIBFg->Draw("CP same");
           TGraphErrors* tgIBF2g = new TGraphErrors(dataNum, gainListData_old, ionBackFlowCorrectedVect1_old, gainErrorListData_old, ionBackFlowCorrectedErrorVect1_old );
           tgIBF2g->SetMarkerStyle(20);
        tgIBF2g->SetMarkerSize(0.3);
           tgIBF2g->SetMarkerColor(3);
           tgIBF2g->Draw("CP same");
        
           TLegend* legendg = new TLegend(0.5,0.65,0.9,0.9);
           //legendg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
           legendg->AddEntry(tgIBF, "Data", "lp");
           legendg->AddEntry(tgIBF2, "Data old", "lp");
           legendg->AddEntry(gr,"Simulation (counting IBF)", "lp");
           //legendg->AddEntry(tgIBFsignal, "Simulation (current signals)", "lp");
           legendg->Draw();


        c5->SaveAs(Form("Figures/IBFCurve-%s-model%d.pdf", gasName.c_str(), modelNum));

     
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    return 0;
}


