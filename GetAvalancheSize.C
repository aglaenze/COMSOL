#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>

#include "Include/Ibf.C"


int GetAvalancheSize() {
    
    bool ionAnalysis = false;
    
    time_t t0 = time(NULL);
    
    LoadStyle();
    
    // input and output files
    const TString path = "rootFiles/Ar-iC4H10/model21/";
    TString fSignalName = path + "signal-250-450-650-900-1010-test.root";
    
    std::cout << fSignalName << std::endl;
    
    
    TCanvas* cv = new TCanvas("cv","cv", 1200, 1000);
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
    
    std::vector<float> *electronStartPointsInput = 0, *electronStartPointsInputX = 0, *electronStartPointsInputY = 0;;
    tAvalanche->SetBranchAddress("electronStartPoints", &electronStartPointsInput);
    tAvalanche->SetBranchAddress("electronStartPointsX", &electronStartPointsInputX);
    tAvalanche->SetBranchAddress("electronStartPointsY", &electronStartPointsInputY);
    
    std::vector<float> electronStartPoints = {}, electronStartPointsX = {}, electronStartPointsY = {};
    std::vector<float> xStart = {}, yStart = {};
    double dmin = -10;
    double dmax = 200;
    TH1F* dElDistribution = new TH1F("dElectrons", "Spread of electrons", int((dmax-dmin)/3.), dmin, dmax);
    const int nAvalanche = tAvalanche->GetEntries();
    for (int k = 0; k < nAvalanche; k++) {
        //for (int k = 4; k < 5; k++) {
        tAvalanche->GetEntry(k);
        electronStartPoints = *electronStartPointsInput;
        electronStartPointsX = *electronStartPointsInputX;
        electronStartPointsY = *electronStartPointsInputY;
        double d = 0;
        int nn = 0;
        for (int j = 0; j< (int)electronStartPoints.size(); j++) {
            if (electronStartPoints[j] > 0.4188 && electronStartPoints[j] < 0.4316) {
                nn++;
                xStart.push_back(electronStartPointsX[j]);
                yStart.push_back(electronStartPointsY[j]);
                //cout << electronStartPointsX[j] << " and " << xStart[-1] << endl;
                //d = 10000*TMath::Sqrt( Square(electronStartPointsX[j]-electronStartPointsX[0]) + Square(electronStartPointsY[j]-electronStartPointsY[0]));
                d = 10000*TMath::Sqrt( Square(xStart[j]-xStart[0]) + Square(yStart[j]-yStart[0]));
                //d = 10000*(xStart[j]-xStart[0]);
                dElDistribution->Fill(d);
            }
        }
        cout << nn << " entries " << endl;
        
        /*
        int yMax = *std::max_element(yStart.begin(), yStart.end());
        int yMin = *std::min_element(yStart.begin(), yStart.end());
        int xMax = *std::max_element(xStart.begin(), xStart.end());
        int xMin = *std::min_element(xStart.begin(), xStart.end());
        double d = 0;
        double x0 = 0, y0 = 0;
        for (int j = 0; j < (int)xStart.size(); j++) {
            x0 = (xMax-xMin)/2.;
            y0 = (yMax-yMin)/2.;
            d = 10000*TMath::Sqrt( Square(xStart[j]-x0) + Square(yStart[j]-y0));
            //d = 10000*(electronStartPointsX[j]-electronStartPointsX[0]);
            dElDistribution->Fill(d);
        }
        xStart.clear();
        yStart.clear();
        */
        electronStartPoints.clear();
        electronStartPointsX.clear();
        electronStartPointsY.clear();
    }
    
    dElDistribution->GetXaxis()->SetTitle("d (um)");
    
    dElDistribution->Draw("hist");
    
    cv->SaveAs("Figures/test-distances.pdf");
    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    return 0;
}


