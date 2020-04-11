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

TF1* GetFitGain(TH1F* h) {
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


int GetGain() {
    
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
    
    double damp = 0.0128, ddrift = 0.5;
    
    std::map <std::string, int> electrode;
    LoadElectrodeMap(modelNum, electrode);
    
    int readoutElectrode = electrode["pad"];

    // Get number of files to look at
    //Int_t num = 1;
    Int_t num = GetNumberOfFiles(path, "convoluted");
    Int_t num2 = int(num/2.);
    if (num/2.> num2) num2+=1;

    TCanvas* c2 = new TCanvas("c2","c2", 600, 300*num2);
    //c2->Divide(2);
    c2->Divide(2, num2);

    Double_t hvMeshList[num];
    Double_t gainList[num], gainList2[num], gainErrorList[num], gainErrorList2[num];
    for (unsigned int k = 0; k < num; ++k) {
        Int_t hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvGemDown = 0, hvGemUp = 0, hvDrift = 0;
        TString fileName;
        if (modelNum == 1) {
            hvMesh = 340+20*k;
            hvDrift = 540+20*k;
            fileName = Form("rootFiles/%s/model%d/fe-spectrum-convoluted-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDrift);
        }
        else if (modelNum >= 2) {
            hvDmDown = 300;
            hvDmUp = 600;
            hvDrift = 800;
            fileName = Form("rootFiles/%s/model%d/fe-spectrum-convoluted-%d-%d-%d.root", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift);
        }
        
        TFile* fGain = TFile::Open(fileName, "READ");
        TH1F* hFeAmplification = (TH1F*)fGain->Get("hFeAmplification");
        TH1F* hFeCharge = (TH1F*)fGain->Get(Form("hFeCharge_%d", readoutElectrode));

        hFeAmplification->Scale(1/hFeAmplification->GetMaximum());
        hFeAmplification->SetMaximum(1.2);
        //hFeAmplification->Rebin(8);
        //hFeAmplification->GetXaxis()->SetRangeUser(2, 10000);
        hFeAmplification->SetLineColor(kBlue + 2);

        TF1* f = GetFitGain(hFeAmplification);
                
        Int_t iBinMax = hFeAmplification->GetMaximumBin();
        Double_t xMax = hFeAmplification->GetXaxis()->GetBinCenter( iBinMax );
        hFeAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
        //hFeAmplification->GetXaxis()->SetRangeUser(0, 40000);
        hvMeshList[k] = hvMesh;
        gainList[k] = f->GetParameter(0);
        gainErrorList[k] = f->GetParError(0);
        
        hFeCharge->Scale(1/hFeCharge->GetMaximum());
        hFeCharge->SetMaximum(1.2);
        hFeCharge->SetLineColor(kBlue + 2);
        TF1* f2 = GetFitGain(hFeCharge);
        hFeCharge->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
        gainList2[k] = f2->GetParameter(0);
        gainErrorList2[k] = f2->GetParError(0);
        
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
        hFeAmplification->Draw("hist");
        f->Draw("same");
        // Add text to frame
        TLatex* txt = new TLatex(0.5*xMax,1,Form("Gain = %.3g #pm %.3f", f->GetParameter(0), f->GetParError(0)));
        TLatex* txt2 = new TLatex(0.5*xMax,0.9,Form("Field ratio = %.2f", (double)hvMesh/(hvDrift-hvMesh)));
        txt->SetTextSize(0.05) ;
        txt2->SetTextSize(0.05) ;
        txt->Draw();
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
        hFeCharge->Draw("hist");
        f2->Draw("same");
        TLatex* txt3 = new TLatex(0.5*xMax,1,Form("Gain = %.3g #pm %.3f", f2->GetParameter(0), f2->GetParError(0)));
        TLatex* txt4 = new TLatex(0.5*xMax,0.9,Form("Field ratio = %.2f", (double)hvMesh/(hvDrift-hvMesh)*(ddrift-damp)/damp));
        txt3->SetTextSize(0.05) ;
        txt4->SetTextSize(0.05) ;
        txt3->Draw();
        txt4->Draw();

    }
    c2->SaveAs(Form("Figures/gains_%s_model%d.pdf", gasName.c_str(), modelNum));

            
    TCanvas* c3 = new TCanvas("c3");
    c3->cd();
    c3->SetGrid();
    c3->SetLogy();
    TGraphErrors* gr = new TGraphErrors(num, hvMeshList, gainList, 0, gainErrorList);
    TGraphErrors* gr2 = new TGraphErrors(num, hvMeshList, gainList2, 0, gainErrorList2);
    gr->SetTitle("Gain curve in the Micromegas");
    gr->GetXaxis()->SetTitle( "V_{mesh}" );
    gr->GetYaxis()->SetTitle( "Gain" );
    gr->GetHistogram()->GetXaxis()->SetLimits(310, 450);
    gr->GetHistogram()->SetMinimum(2.e2);   // along Y axis
    gr->GetHistogram()->SetMaximum(8.e4);   // along Y axis
    gr->SetMarkerStyle(20);
    gr->Draw("ACP");
    gr->Draw("ACP same");
    
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
    Double_t hvMeshListData[dataNum], hvDriftListData[dataNum], gainListData[dataNum], gainErrorListData_old[dataNum], gainListData_old[dataNum], gainErrorListData[dataNum];
    LoadGainData(gasName, dataNum, hvMeshListData, hvDriftListData, gainListData, gainErrorListData, gainListData_old, gainErrorListData_old);
        
    TGraphErrors* grd = new TGraphErrors(dataNum, hvMeshListData, gainListData, 0, gainErrorListData);
    grd->SetMarkerStyle(20);
    grd->Draw("CP same");
    TGraphErrors* grd2 = new TGraphErrors(dataNum, hvMeshListData, gainListData_old, 0, gainErrorListData_old);
    grd2->SetMarkerStyle(20);
    grd2->Draw("CP same");
    
    TF1* fExp2 = new TF1( "FitFunctionExp2", FitFunctionExp, hvMeshListData[0]-5, hvMeshListData[dataNum-1]+5, 2 );
    fExp2->SetParName( 0, "const" );
    fExp2->SetParName( 1, "slope" );
    fExp2->SetLineColor(3);
    grd->Fit( fExp2, "0" );
    fExp2->Draw("same");
    PutText( 0.2, 0.8, Form( "y_{data} = exp( %.3g + %.3g x)", fExp2->GetParameter(0), fExp2->GetParameter(1) ) );
        
    TF1* fExp3 = new TF1( "FitFunctionExp3", FitFunctionExp, hvMeshListData[0]-5, hvMeshListData[dataNum-1]+5, 2 );
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

     
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    return 0;
}


