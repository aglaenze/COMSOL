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
    const bool drawIbf = false;
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
        Double_t hvMeshList[num];
        Double_t gainList[num];
        Double_t gainErrorList[num];
        for (unsigned int k = 0; k < num; ++k) {
            Int_t Vmesh = 320 + k*20;
            if (modelNum == 4) Vmesh = 350 + k*10;
            TString fileName = Form("rootFiles/%s/model%d/Fe_spectrum_convoluted_%dV.root", gasName.c_str(), modelNum, Vmesh);
            TFile* fGain = TFile::Open(fileName, "READ");
            TH1F* hGain = (TH1F*)fGain->Get("hFeElectrons");
            hGain->SetTitle(Form("Number of secondary electrons for V_{mesh} = %d V with a simulated Fe source", Vmesh));
            //hGain->Rebin(8);
            //hGain->GetXaxis()->SetRangeUser(2, 10000);
            //hGain->SetFillColor(kBlue + 2);
            hGain->SetLineColor(kBlue + 2);
            
            Int_t iBinMax = hGain->GetMaximumBin();
            Int_t xMax = hGain->GetXaxis()->GetBinCenter( iBinMax );
            
            /*
            std::cout << "Maximum = " << hGain->GetMaximum() << std::endl;
            std::cout << "RMS = " << hGain->GetRMS() << std::endl;
            std::cout << "Mean = " << hGain->GetMean() << std::endl;
            std::cout << "xMax = " << xMax << std::endl;
            std::cout << "iBinMax = " << iBinMax << std::endl;
            */
            //std::cout << "\n\nNormalisation = " << hGain->GetMaximum()*hGain->GetRMS() << std::endl;
            
            Int_t fitRangeMin = xMax - 1.0 * hGain->GetRMS();
            Int_t fitRangeMax = xMax + 1.0 * hGain->GetRMS();

            /*
            TF1* f = new TF1( "FitFunction", FitFunctionCrystalBall, fitRangeMin, fitRangeMax, 5);
            f->SetParNames("normalisation", "alpha", "n", "Sigma", "Mean");
            f->SetParameters(hGain->GetMaximum()*hGain->GetRMS(), 3, 3, hGain->GetRMS(), xMax);
             */
            
            TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
            f->SetParNames("Mean", "Sigma", "Amplitude");
            f->SetParameters(xMax, hGain->GetRMS(), hGain->GetMaximum());
            
            hGain->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
            
            
            c2->cd(k+1);
            //c2->SetLogy();
            hGain->GetXaxis()->SetRangeUser(0, xMax + 3*hGain->GetRMS());
            hGain->Draw();
            f->Draw("same");
            hvMeshList[k] = Vmesh;
            gainList[k] = f->GetParameter(0);
            gainErrorList[k] = f->GetParError(0);
            //std::cout << "Gain = " << gainList[k] << std::endl;
        }
        c2->SaveAs(Form("Figures/gains_%s_model%d.pdf", gasName.c_str(), modelNum));

            TCanvas* c3 = new TCanvas("c3");
            c3->cd();
            c3->SetGrid();
            c3->SetLogy();
            TGraphErrors* gr = new TGraphErrors(num, hvMeshList, gainList, 0, gainErrorList);
            gr->SetTitle("Gain curve in the Micromegas");
            gr->GetXaxis()->SetTitle( "V_{mesh}" );
            gr->GetYaxis()->SetTitle( "Gain" );
            gr->GetHistogram()->GetXaxis()->SetLimits(310, 450);
            gr->GetHistogram()->SetMinimum(2.e2);   // along Y axis
            gr->GetHistogram()->SetMaximum(8.e4);   // along Y axis
            gr->SetMarkerStyle(20);
            gr->Draw("ACP");

            // Fit the gain curve with an exponential function
            TF1* fExp = new TF1( "FitFunctionExp", FitFunctionExp, hvMeshList[0]-5, hvMeshList[num-1]+5, 2 );
            fExp->SetParName( 0, "const" );
            fExp->SetParName( 1, "slope" );
        
            gr->Fit( fExp, "0");
            fExp->SetLineColor(2);
            fExp->Draw("same");
            PutText( 0.2, 0.7, Form( "y_{simulation} = exp( %.3g + %.3g x)", fExp->GetParameter(0), fExp->GetParameter(1) ) );
        
        // Same with data
        const Int_t dataNum = dataQuantity(gasName);
        Double_t gainListData[dataNum], gainErrorListData[dataNum], hvMeshListData[dataNum];
        LoadGainData(gasName, dataNum, hvMeshListData, gainListData, gainErrorListData);
        
            TGraphErrors* grd = new TGraphErrors(dataNum, hvMeshListData, gainListData, 0, gainErrorListData);
            grd->SetMarkerStyle(20);
            grd->Draw("CP same");
        
            TF1* fExp2 = new TF1( "FitFunctionExp2", FitFunctionExp, hvMeshListData[0]-5, hvMeshListData[dataNum-1]+5, 2 );
            fExp2->SetParName( 0, "const" );
            fExp2->SetParName( 1, "slope" );
            fExp2->SetLineColor(3);
            grd->Fit( fExp2, "0" );
            fExp2->Draw("same");
            PutText( 0.2, 0.8, Form( "y_{data} = exp( %.3g + %.3g x)", fExp2->GetParameter(0), fExp2->GetParameter(1) ) );
        
            TLegend* legend = new TLegend(0.7,0.2,0.9,0.4);
            //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
            legend->AddEntry(fExp2,"Data","lp");
            legend->AddEntry(fExp,"Simulation","lp");
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
        Double_t hvMeshList[num];
        Double_t ibfList[num];
        Double_t ibfErrorList[num];
        
        for (unsigned int k = 0; k < num; ++k) {
        //Int_t k = 3;
            Int_t Vmesh = 320 + k*20;
            TString fileName = path + Form("ibf_%dV.root", Vmesh);
            TFile* fIbf = TFile::Open(fileName, "READ");
            TH1F* hIbf = (TH1F*)fIbf->Get("hIbf");
            hIbf->SetTitle(Form("IBF for for V_{mesh} = %d V", Vmesh));
            TF1* f = new TF1( "FitFunction", FitGauss, 0, 10, 3);
            f->SetParNames("Mean", "Sigma", "Amplitude");
            f->SetParameters(hIbf->GetMean(), hIbf->GetRMS(), hIbf->GetMaximum());
            hIbf->Fit(f, "0");
            hIbf->Rebin(2);
            c4->cd(k+1);
        //c4->cd();
            hIbf->Draw();
            f->Draw("same");
            hvMeshList[k] = Vmesh;
            ibfList[k] = f->GetParameter(0);
            ibfErrorList[k] = f->GetParError(0);
            
        }
        c4->SaveAs(Form("Figures/ibfResults_%s_model%d.pdf", gasName.c_str(), modelNum));
        
        Int_t numS = GetNumberOfFiles(path, "signal");
        Double_t hvMeshListS[numS];
        Double_t ibfSignalList[numS];

        for (unsigned int k = 0; k < numS; ++k) {
        //Int_t k = 3;
            Int_t Vmesh = 320 + k*20;
            TString fileName = path + Form("signal_%dV.root", Vmesh);
            TFile* fSignal = TFile::Open(fileName, "READ");
            TH1F* hInt = (TH1F*)fSignal->Get("hIntMesh");
            TH1F* hIntd = (TH1F*)fSignal->Get("hIntDrift");
            hvMeshListS[k] = Vmesh;
            Int_t last = hInt->GetNbinsX();
            ibfSignalList[k] =  100.*hIntd->GetBinContent(last)/hInt->GetBinContent(last);
        }
        

        TCanvas* c5 = new TCanvas("c5");
        c5->cd();
        c5->SetGrid();
        TGraphErrors* gr = new TGraphErrors(num, hvMeshList, ibfList, 0, ibfErrorList);
        gr->SetTitle("IBF curve in the Micromegas");
        gr->GetXaxis()->SetTitle( "V_{mesh}" );
        gr->GetYaxis()->SetTitle( "IBF (%)" );
        //gr->GetXaxis()->SetLimits(0, 4);  // along X axis
        gr->GetHistogram()->SetMinimum(0.);   // along Y axis
        gr->GetHistogram()->SetMaximum(4.);   // along Y axis
        gr->SetMarkerStyle(20);
        gr->GetXaxis()->SetLimits(315, 445);
        gr->Draw("ACP");
        
        TGraph* tgIBFsignal = new TGraph(numS, hvMeshListS, ibfSignalList );
        tgIBFsignal->SetMarkerStyle(20);
        tgIBFsignal->SetMarkerColor(3);
        tgIBFsignal->Draw("CP same");
        
     // Same with data
     const Int_t dataNum = dataQuantity(gasName);
     Double_t hvMeshListIBF1[dataNum], ionBackFlowCorrectedVect1[dataNum], ionBackFlowCorrectedErrorVect1[dataNum];
     LoadGainData(gasName, dataNum, hvMeshListIBF1, ionBackFlowCorrectedVect1, ionBackFlowCorrectedErrorVect1);
        
        TGraphErrors* tgIBF = new TGraphErrors(dataNum, hvMeshListIBF1, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
        tgIBF->SetMarkerStyle(20);
        tgIBF->SetMarkerColor(2);
        tgIBF->Draw("CP same");
        
        TLegend* legend = new TLegend(0.5,0.65,0.9,0.9);
        //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
        legend->AddEntry(tgIBF,"Data","lp");
        legend->AddEntry(gr,"Simulation (counting IBF)","lp");
        legend->AddEntry(tgIBFsignal, "Simulation (current signals)", "lp");
        legend->Draw();

        c5->SaveAs(Form("Figures/IBFCurve_%s_model%d.pdf", gasName.c_str(), modelNum));
    }
    
    
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    return 0;
}


