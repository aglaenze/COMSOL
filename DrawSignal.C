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



int DrawSignal() {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    //____________________

    time_t t0 = time(NULL);
    gStyle->SetOptStat(0);
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
     
    //if (gainCurve) {
    // Get number of files to look at
    //Int_t num = 1;
    Int_t num = GetNumberOfFiles(path, "signal");
    /*Int_t num2 = int(num/2.);
    if (num/2.> num2) num2+=1;
    TCanvas* c2 = new TCanvas("c2");
    c2->Divide(2, num2);
    */
    Int_t electrodeNum = 0;
    if (modelNum == 1) electrodeNum = 3;
    else if (modelNum >= 2 && modelNum < 5) electrodeNum = 4;
    else if (modelNum >= 5 && modelNum < 10) electrodeNum = 5;
    TCanvas* cvf = new TCanvas("cvf");
    cvf->Divide(electrodeNum, num);
    
    Double_t hvMmList[num];
    Double_t gainList[num];
    Double_t gainErrorList[num];
    for (unsigned int k = 0; k < num; ++k) {
        Int_t hvMm = 0, hvMmDown = 0, hvMmUp = 0, hvDrift = 0;
        TString fileName;
        if (modelNum == 1) {
            hvMm = 340+20*k;
            hvDrift = 540+20*k;
            fileName = path + Form("signal-%d-%d.root", hvMm, hvDrift);
        }
        else if (modelNum > 6 && modelNum < 10) {
            hvMmDown = 300;
            hvMmUp = 600;
            hvDrift = 800;
            fileName = path + Form("signal-%d-%d-%d.root", hvMmDown, hvMmUp, hvDrift);
        }
        else if (modelNum >= 10 && modelNum < 13) {
            hvMm = 300;
            hvMmDown = 400;
            hvMmUp = 700;
            hvDrift = 900;
            fileName = path + Form("signal-%d-%d-%d-%d.root", hvMm, hvMmDown, hvMmUp, hvDrift);
        }
        // import data
        TFile* fSignal = TFile::Open(fileName, "READ");
  
        for (int j = 0; j < electrodeNum; j++) {
        //for (int j = 0; j < 1; j++) {
            TTree* tSignal = (TTree*)fSignal->Get(Form("tSignal_V%d", j+1));
            /*
            int nEntries = tSignal->GetEntries();
            double time[nEntries], ct[nEntries], ce[nEntries], ci[nEntries];
            // following parameters are taken from signal.C
            const double tStep = 0.1;  //ns
            const double tStart = 0.;
            const double tEnd = tSignal->GetMaximum("time");
            std::cout << "tEnd = " << tEnd << std::endl;
            const int nTimeBins = (tEnd - tStart)/tStep;
        
            const double iMin =  tSignal->GetMinimum("electronCurrent");
            const double iMax =  tSignal->GetMaximum("electronCurrent");
            std::cout << "iMax = " << iMax << std::endl;
            const double iStep = 1.e-2;
            const int nCurrentBins = abs((iMax-iMin)/iStep);
            std::cout << "nCurrentBins = " << nCurrentBins << std::endl;
            TH2F* hESignal = new TH2F("hESignal", "Signal = f(t)", nTimeBins, tStart, tEnd, nCurrentBins, iMin, iMax);
            //tSignal->Draw("electronCurrent:time>>hESignal");
             */
             
                /*
                TH1F* hESignal = new TH1F("hESignal", "Signal = f(t)", nTimeBins, tStart, tEnd);
                double* rowContent;
                //TH1F hESignal;
                //tSignal->Draw("electronCurrent:time","","hist");
                for (int l = 0; l < nEntries; l++) {
                    tSignal->GetEntry(irow);
                    rowContent = tSignal->GetArgs();
                    hESignal->Fill(rowContent[0], rowContent[1]);
                }
                 */
              /*
            TCanvas* cv = new TCanvas(Form("cv%d", j+1));
            hESignal->SetLineColor(kRed);
            hESignal->SetMarkerStyle(kPlus);
            hESignal->SetLineWidth(10);
            //hESignal->Draw("HIST C");
               */
            //hESignal->Draw("box");
            //cv->SaveAs(Form("Figures/SignalTestV%d.pdf", j+1));
            //return;
            
            cvf->cd(j*num+k+1);
            std::cout << j*num+k+1 << std::endl;
            // find how to change colors
            tSignal->Draw("electronCurrent");
            tSignal->Draw("ionCurrent","", "same");
            tSignal->Draw("totalCurrent","", "same");
        }
     }
     
    cvf->SaveAs("Figures/ElectronCurrents.pdf");
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    return 0;
}


