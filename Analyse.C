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



Double_t FitFunctionCrystalBall( Double_t* x, Double_t* par ) { //(x, alpha, n sigma, mu)
    return par[0]*ROOT::Math::crystalball_pdf( x[0], par[1], par[2], par[3], par[4] ); }

Double_t FitGauss( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
    return  par[2]*TMath::Gaus( x[0], par[0], par[1]); }

//____________________________________________
Double_t FitFunctionExp( Double_t* x, Double_t* par ) {
    return TMath::Exp( par[0] + par[1]*x[0] ); }


int Analyse() {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    const bool gainCurve = true;
    const bool drawIbf = true;
    //____________________

    time_t t0 = time(NULL);
    gStyle->SetOptStat(0);
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
     
    if (gainCurve) {
        // Get number of files to look at
        //Int_t num = 3;
        Int_t num = GetNumberOfFiles(path, "convoluted");
        Int_t num2 = int(num/2.);
        if (num/2.> num2) num2+=1;
        TCanvas* c2 = new TCanvas("c2");
        c2->Divide(2, num2);
        Double_t hvMmList[num];
        Double_t gainList[num];
        Double_t gainErrorList[num];
        for (unsigned int k = 0; k < num; ++k) {
            Int_t hvMm = 0, hvMmDown = 0, hvMmUp = 0, hvDrift = 0;
            TString fileName;
            if (modelNum == 1) {
                hvMm = 340+20*k;
                hvDrift = 540+20*k;
                fileName = Form("rootFiles/%s/model%d/Fe-spectrum-convoluted-%d-%d.root", gasName.c_str(), modelNum, hvMm, hvDrift);
            }
            else if (modelNum > 6) {
                hvMmDown = 300;
                hvMmUp = 600;
                hvDrift = 800;
                fileName = Form("rootFiles/%s/model%d/Fe-spectrum-convoluted-%d-%d-%d.root", gasName.c_str(), modelNum, hvMmDown, hvMmUp, hvDrift);
            }
            TFile* fGain = TFile::Open(fileName, "READ");
            TH1F* hGain = (TH1F*)fGain->Get("hFeElectrons");
            /*
            int nBins = hGain1->GetNbinsX();
            Double_t xMaxh = hGain1->GetXaxis()->GetBinCenter( nBins);
            TH1F* hGain = new TH1F("", "", 40000*nBins/xMaxh, 0, 40000 );
            //std::cout << hGain1->GetEntries() << std::endl;
            for (int i=0; i<hGain1->GetEntries(); i++) {
                Double_t value = hGain1->GetBinContent( i );
                Double_t x = hGain1->GetXaxis()->GetBinCenter( i);
                hGain->Fill(x, value);
            }
             */
            hGain->SetTitle(Form("Number of secondary electrons for V_{mesh} = %d V with a simulated Fe source", hvMm));
            hGain->Scale(1/hGain->GetMaximum());
            hGain->SetMaximum(1.2);
            //hGain->Rebin(8);
            //hGain->GetXaxis()->SetRangeUser(2, 10000);
            hGain->SetLineColor(kBlue + 2);
            
            Int_t iBinMax = hGain->GetMaximumBin();
            Int_t xMax = hGain->GetXaxis()->GetBinCenter( iBinMax );
            
            Int_t fitRangeMin = xMax - 0.6 * hGain->GetRMS();
            Int_t fitRangeMax = xMax + 0.6 * hGain->GetRMS();

            TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
            f->SetParNames("Mean", "Sigma", "Amplitude");
            f->SetParameters(xMax, hGain->GetRMS(), hGain->GetMaximum());
            
            hGain->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
                        
            c2->cd(k+1);
            hGain->GetXaxis()->SetRangeUser(0, xMax + 3*hGain->GetRMS());
            //hGain->GetXaxis()->SetRangeUser(0, 40000);
            hGain->Draw("hist");
            f->Draw("same");
            hvMmList[k] = hvMm;
            gainList[k] = f->GetParameter(0);
            gainErrorList[k] = f->GetParError(0);
            //std::cout << "Gain = " << gainList[k] << std::endl;
            //c2->SaveAs(Form("Figures/gains_%s_model%d_%d.pdf", gasName.c_str(), modelNum, hvMm));
        }
        c2->SaveAs(Form("Figures/gains_%s_model%d.pdf", gasName.c_str(), modelNum));

            
        TCanvas* c3 = new TCanvas("c3");
        c3->cd();
        c3->SetGrid();
        c3->SetLogy();
        TGraphErrors* gr = new TGraphErrors(num, hvMmList, gainList, 0, gainErrorList);
        gr->SetTitle("Gain curve in the Micromegas");
        gr->GetXaxis()->SetTitle( "V_{mesh}" );
        gr->GetYaxis()->SetTitle( "Gain" );
        gr->GetHistogram()->GetXaxis()->SetLimits(310, 450);
        gr->GetHistogram()->SetMinimum(2.e2);   // along Y axis
        gr->GetHistogram()->SetMaximum(8.e4);   // along Y axis
        gr->SetMarkerStyle(20);
        gr->Draw("ACP");

        // Fit the gain curve with an exponential function
        TF1* fExp = new TF1( "FitFunctionExp", FitFunctionExp, hvMmList[0]-5, hvMmList[num-1]+5, 2 );
        fExp->SetParName( 0, "const" );
        fExp->SetParName( 1, "slope" );
        gr->Fit( fExp, "0");
        fExp->SetLineColor(2);
        fExp->Draw("same");
        PutText( 0.2, 0.7, Form( "y_{simulation} = exp( %.3g + %.3g x)", fExp->GetParameter(0), fExp->GetParameter(1) ) );
        
        // Same with data
        const Int_t dataNum = dataQuantity(gasName);
        Double_t hvMmListData[dataNum], hvDriftListData[dataNum], gainListData[dataNum], gainErrorListData_old[dataNum], gainListData_old[dataNum], gainErrorListData[dataNum];
        LoadGainData(gasName, dataNum, hvMmListData, hvDriftListData, gainListData, gainErrorListData, gainListData_old, gainErrorListData_old);
        
        TGraphErrors* grd = new TGraphErrors(dataNum, hvMmListData, gainListData, 0, gainErrorListData);
        grd->SetMarkerStyle(20);
        grd->Draw("CP same");
        TGraphErrors* grd2 = new TGraphErrors(dataNum, hvMmListData, gainListData_old, 0, gainErrorListData_old);
        grd2->SetMarkerStyle(20);
        grd2->Draw("CP same");
    
        TF1* fExp2 = new TF1( "FitFunctionExp2", FitFunctionExp, hvMmListData[0]-5, hvMmListData[dataNum-1]+5, 2 );
        fExp2->SetParName( 0, "const" );
        fExp2->SetParName( 1, "slope" );
        fExp2->SetLineColor(3);
        grd->Fit( fExp2, "0" );
        fExp2->Draw("same");
        PutText( 0.2, 0.8, Form( "y_{data} = exp( %.3g + %.3g x)", fExp2->GetParameter(0), fExp2->GetParameter(1) ) );
        
        TF1* fExp3 = new TF1( "FitFunctionExp3", FitFunctionExp, hvMmListData[0]-5, hvMmListData[dataNum-1]+5, 2 );
        fExp3->SetParName( 0, "const" );
        fExp3->SetParName( 1, "slope" );
        fExp3->SetLineColor(4);
        grd2->Fit( fExp3, "0" );
        fExp3->Draw("same");
        PutText( 0.2, 0.75, Form( "y_{data old} = exp( %.3g + %.3g x)", fExp3->GetParameter(0), fExp3->GetParameter(1) ) );
    
        TLegend* legend = new TLegend(0.7,0.2,0.9,0.4);
        //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend->AddEntry(fExp2,"Data", "l");
        legend->AddEntry(fExp3,"Data (old)", "l");
        legend->AddEntry(fExp,"Simulation", "l");
        legend->Draw();
        
        c3->SaveAs(Form("Figures/GainCurve_%s_model%d.pdf", gasName.c_str(), modelNum));
    }
    
    
    if (drawIbf) {
        
        // Get number of files to look at
        //Int_t num = 3;
        Int_t num = GetNumberOfFiles(path, "ibf");
        Int_t num2 = int(num/2.);
        if (num/2.> num2) num2+=1;
        TCanvas* c4 = new TCanvas("c4");
        c4->Divide(2, num2);
        Double_t hvMmList[num], hvDriftList[num], hvRatioList[num];
        Double_t ibfList[num], ibfErrorList[num];
        
        double damp = 0.0128, ddrift = 0.5;
        
        for (unsigned int k = 0; k < num; ++k) {
        //Int_t k = 3;
            Int_t hvMm = 0, hvMmDown = 0, hvMmUp = 0, hvDrift = 0;
            TString fileName;
            if (modelNum == 1) {
                hvMm = 340+20*k;
                hvDrift = 540+20*k;
                fileName = path + Form("ibf-%d-%d.root", hvMm, hvDrift);
            }
            else if (modelNum > 6) {
                hvMmDown = 300;
                hvMmUp = 600;
                hvDrift = 800;
                fileName = path + Form("ibf-%d-%d-%d.root", hvMmDown, hvMmUp, hvDrift);
            }
            TFile* fIbf = TFile::Open(fileName, "READ");
            TH1F* hIbf = (TH1F*)fIbf->Get("hIbf");
            hIbf->SetTitle(Form("IBF for for V_{mesh} = %d V", hvMm));
            TF1* f = new TF1( "FitFunction", FitGauss, 0, 10, 3);
            f->SetParNames("Mean", "Sigma", "Amplitude");
            f->SetParameters(hIbf->GetMean(), hIbf->GetRMS(), hIbf->GetMaximum());
            hIbf->Fit(f, "0");
            //hIbf->Rebin(2);
            c4->cd(k+1);
            hIbf->Draw();
            f->Draw("same");
            hvMmList[k] = hvMm;
            hvDriftList[k] = hvDrift;
            hvRatioList[k] = hvMm/(hvDrift-hvMm) * (ddrift-damp)/damp;
            ibfList[k] = f->GetParameter(0);
            ibfErrorList[k] = f->GetParError(0);
            std::cout << hvMm << " " << hvDrift << " " << hvRatioList[k] << std::endl;
            
        }
        return 0;
        c4->SaveAs(Form("Figures/ibfResults-%s-model%d.pdf", gasName.c_str(), modelNum));
        
        Int_t numS = GetNumberOfFiles(path, "signal");
        Double_t hvMmListS[numS];
        Double_t ibfSignalList[numS];

        for (unsigned int k = 0; k < numS; ++k) {
        //Int_t k = 3;
            Int_t hvMmS = 320 + k*20;
            TString fileName = path + Form("signal_%dV.root", hvMmS);
            TFile* fSignal = TFile::Open(fileName, "READ");
            TH1F* hInt = (TH1F*)fSignal->Get("hIntMesh");
            TH1F* hIntd = (TH1F*)fSignal->Get("hIntDrift");
            hvMmListS[k] = hvMmS;
            Int_t last = hInt->GetNbinsX();
            ibfSignalList[k] =  100.*hIntd->GetBinContent(last)/hInt->GetBinContent(last);
        }
        

        TCanvas* c5 = new TCanvas("c5", "c5", 200, 500);
        c5->SetGrid();
        c5->Divide(1,3);
        
        c5->cd(1);
        TGraphErrors* gr = new TGraphErrors(num, hvMmList, ibfList, 0, ibfErrorList);
        gr->SetTitle("IBF curve in the Micromegas");
        gr->GetXaxis()->SetTitle( "V_{mesh}" );
        gr->GetYaxis()->SetTitle( "IBF (%)" );
        //gr->GetXaxis()->SetLimits(0, 4);  // along X axis
        gr->GetHistogram()->SetMinimum(0.);   // along Y axis
        gr->GetHistogram()->SetMaximum(4.);   // along Y axis
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.3);
        //gr->GetXaxis()->SetLimits(hvMmList[0]-5, hvMmList[num]+5);
        gr->Draw("ACP");
        
        TGraph* tgIBFsignal = new TGraph(numS, hvMmListS, ibfSignalList );
        tgIBFsignal->SetMarkerStyle(20);
        tgIBFsignal->SetMarkerSize(0.3);
        tgIBFsignal->SetMarkerColor(3);
        //tgIBFsignal->Draw("CP same");
        
        // Same with data
        const Int_t dataNum = dataQuantity(gasName);
        Double_t hvMmListIBF1[dataNum], hvDriftListIBF1[dataNum], ionBackFlowCorrectedVect1[dataNum], ionBackFlowCorrectedErrorVect1[dataNum], ionBackFlowCorrectedVect1_old[dataNum], ionBackFlowCorrectedErrorVect1_old[dataNum];
        LoadIbfData(gasName, dataNum, hvMmListIBF1, hvDriftListIBF1, ionBackFlowCorrectedVect1, ionBackFlowCorrectedErrorVect1, ionBackFlowCorrectedVect1_old, ionBackFlowCorrectedErrorVect1_old);
        
        Double_t ratioListData[dataNum];
        for (int i = 0; i < dataNum; i++) {ratioListData[i] = hvMmListIBF1[i]/(hvDriftListIBF1[i]-hvMmListIBF1[i]) * (ddrift-damp)/damp;}
        
        PrintList("ratioListData", ratioListData, num);
        PrintList("hvRatioList", hvRatioList, num);
        return 0;
         
        TGraphErrors* tgIBF = new TGraphErrors(dataNum, hvMmListIBF1, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
        tgIBF->SetMarkerStyle(20);
        tgIBF->SetMarkerSize(0.3);
        tgIBF->SetMarkerColor(2);
        tgIBF->Draw("CP same");
        TGraphErrors* tgIBF2 = new TGraphErrors(dataNum, hvMmListIBF1, ionBackFlowCorrectedVect1_old, 0, ionBackFlowCorrectedErrorVect1_old );
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
           gr2->SetTitle("IBF curve in the Micromegas");
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


        c5->SaveAs(Form("Figures/IBFCurve-%s-model%d.pdf", gasName.c_str(), modelNum));
     }
     
     
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    return 0;
}


