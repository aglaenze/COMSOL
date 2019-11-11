#include <iostream>
#include <fstream>
#include <cmath>
#include <dirent.h>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "_Utils.C"
#include "parameters.C"

#include <cmath>

/* Creates histograms that are the convolution of Fe55 spectrum histogram (number of primaries created by 1 photon) and secondaries histogram (number of secondaries created by one primary in the drift region, after amplification)
 */


int Convolute() {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    //____________________
    
    
    time_t t0 = time(NULL);
    
    gStyle->SetOptStat(0);

    TFile fFe(Form("rootFiles/%s/spectrum_Fe55.root", gasName.c_str()));
    TH1F* hFe = (TH1F*)fFe.Get("hElectrons");
     
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    Int_t num = GetNumberOfFiles(path, "gain");
    
    for (unsigned int k = 0; k < num; ++k) {
        Int_t Vmesh = 320 + k*20;

        TString fGainName = path+Form("gain_%dV", Vmesh);
        TFile fGain(fGainName+".root");
        TH1F* hGain = (TH1F*)fGain.Get("hElectrons");
        
        Int_t nGain = hGain->GetEntries();
        Int_t nFe = hFe->GetEntries();
        //Int_t nFe = 228; // nPrimary for Ar-iC_{4}H_{10} 95/5
        std::cout << "Number of entries in hGain = " << nGain << std::endl;
        std::cout << "Number of entries in hFe = " << nFe << std::endl;
        
        Int_t nPrimaryTh = GetPrimary(gasName);

        const Int_t nBins = int(nGain/6);
        TH1F* hFeElectrons = new TH1F("hFeElectrons", "Number of secondary electrons with Fe source", nBins, 0, hGain->GetMaximumBin()*40);
        
        for (unsigned int i = 0; i < 100000; ++i) {
            Int_t nPrim = hFe->GetRandom();
            //std::cout << "Nprim = " << nPrim << std::endl;
            Int_t gtot = 0;
            for (unsigned int j = 0; j < nPrim; ++j) {
                Int_t gain = hGain->GetRandom();
                //std::cout << "gain = " << gain << std::endl;
                //std::cout << "gtot = " << gtot << std::endl;
                gtot += gain;
                //hFeElectrons->Fill(gain);
            }
            hFeElectrons->Fill(gtot/nPrimaryTh);
            //std::cout << gtot << std::endl;
        }
        
        TCanvas* c1 = new TCanvas();
        c1->cd();
        
        hFeElectrons->SetXTitle("# secondary electrons");
        hFeElectrons->SetYTitle("# counts");
        hFeElectrons->Draw();
        //c1->SaveAs(Form("Convolution_%s.pdf", gasName.c_str()));
        
        // Write convolution histograms in root files
        TFile* f = new TFile(path+Form("Fe_spectrum_convoluted_%dV.root", Vmesh), "RECREATE");
        hFeElectrons->Write();
        f->Close();
    }

    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    return 0;
}


