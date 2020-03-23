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

#include <cmath>

/* Creates histograms that are the convolution of Fe55 spectrum histogram (number of primaries created by 1 photon) and secondaries histogram (number of secondaries created by one primary in the drift region, after amplification)
 */


int Convolute() {
    
    //______________________
    // variables
    //std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    std::string gasName = "Ar-iC4H10";
    const int modelNum = 1;
    //____________________
    
    
    time_t t0 = time(NULL);
    
    gStyle->SetOptStat(0);

    TFile fFe(Form("rootFiles/%s/spectrum_Fe55.root", gasName.c_str()));
    TH1F* hFe = (TH1F*)fFe.Get("hElectrons");
     
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    Int_t num = GetNumberOfFiles(path, "signal");
    num = 1;
    
    Int_t nPrimaryTh = GetPrimary(gasName);
    //Int_t nPrimaryTh = 600;
    
    TString fSignalName, fChargeName, fOutputName;
    
    Int_t hvMm = 0, hvDmDown = 0, hvDmUp = 0, hvDrift = 0;
    hvMm = 340;
    hvDrift = 540;
    fSignalName = path+Form("signal-%d-%d.root", hvMm, hvDrift);
    fChargeName = path+Form("charges-%d-%d.root", hvMm, hvDrift);
    fOutputName = path+Form("Fe-spectrum-convoluted-%d-%d.root", hvMm, hvDrift);
    
    std::map <std::string, int> electrode;
    if (modelNum==1) {
        electrode["mesh"] = 2;
        electrode["drift"] = 3;
        electrode["pad"] = 4;
    }
    
    //for (unsigned int k = 0; k < num; ++k) {
    /*
        Int_t hvMm = 0, hvDmDown = 0, hvDmUp = 0, hvDrift = 0;
        if (modelNum == 1) {
            hvMm = 340+20*k;
            hvDrift = 540+20*k;
            fsignalName = path+Form("signal-%d-%d", hvMm, hvDrift);
            fOutputName = path+Form("Fe-spectrum-convoluted-%d-%d.root", hvMm, hvDrift);
        }
        else if (modelNum >= 2 && modelNum < 5) {
            hvDmDown = 300;
            hvDmUp = 600;
            hvDrift = 800;
            fsignalName = path+Form("signal-%d-%d-%d", hvDmDown, hvDmUp, hvDrift);
            fOutputName = path+Form("Fe-spectrum-convoluted-%d-%d-%d.root", hvDmDown, hvDmUp, hvDrift);
        }
        else if (modelNum == 5) {
            hvMm = 300;
            hvDmDown = 400;
            hvDmUp = 700;
            hvDrift = 900;
            fsignalName = path+Form("signal-%d-%d-%d-%d", hvMm, hvDmDown, hvDmUp, hvDrift);
            fOutputName = path+Form("Fe-spectrum-convoluted-%d-%d-%d-%d.root", hvMm, hvDmDown, hvDmUp, hvDrift);
        }

        TFile fsignal(fSignalName + ".root");
        TH1F* hsignal = (TH1F*) fsignal.Get("hElectrons");
        //hsignal->Rebin(15);
     */
    
    int readoutElectrode = electrode["pad"];
    
    TFile fSignal(fSignalName);
    TTree* tGain = (TTree*) fSignal.Get("tGain");
    Int_t gain;
    tGain->SetBranchAddress("secondaryElectrons", &gain);
    
    TFile fCharge(fChargeName);
    TTree* tCharge = (TTree*) fCharge.Get(Form("tCharges_%d", readoutElectrode));
    Double_t nTotal;
    tCharge->SetBranchAddress("nTotal", &nTotal);
        
        Int_t nGain = tGain->GetEntries();
        Int_t nCharge = tCharge->GetEntries();
        Int_t nFe = hFe->GetEntries();
        std::cout << "Number of entries in tGain = " << nGain << std::endl;
        std::cout << "Number of entries in hFe = " << nFe << std::endl;
    std::cout << "Number of entries in tCharge = " << nCharge << std::endl;

        //const Int_t nBins = int(nGain/4);
        const Int_t nBins = int(tGain->GetMaximum("secondaryElectrons")/4);
    const Int_t nBins2 = int(tCharge->GetMaximum("nTotal")/4);
    
        TH1F* hFeElectrons = new TH1F("hFeElectrons", "Number of secondary electrons with Fe source", nBins, 0, nBins*4 );
    TH1F* hFeSignal = new TH1F("hFeSignal", "Signal with Fe source", nBins2, 0, nBins2*4 );
        
        //std::cout << "maximum bin = " << hsignal->GetMaximumBin() << std::endl;
        
        for (unsigned int i = 0; i < 100000; ++i) {
            Int_t nPrim = hFe->GetRandom();
            //std::cout << "\nNprim = " << nPrim << std::endl;
            Double_t gtot = 0, gtot2 = 0;
            for (unsigned int j = 0; j < nPrim; ++j) {
                int r = rand() % nGain;
                tGain->GetEntry(r);
                gtot += gain;
                int r2 = rand() % nCharge;
                //return 0;
                tCharge->GetEntry(r2);
                gtot2 += nTotal;
            }
            hFeElectrons->Fill(gtot/nPrimaryTh);
            hFeSignal->Fill(gtot2/nPrimaryTh);
            //std::cout << gtot << std::endl;
        }
        
        TCanvas* c1 = new TCanvas();
        c1->cd();
        
        hFeElectrons->SetXTitle("# secondary electrons");
        hFeElectrons->SetYTitle("# counts");
        hFeElectrons->Draw();
        //c1->SaveAs(Form("Convolution_%s.pdf", gasName.c_str()));
    
    hFeSignal->SetXTitle("# secondary electrons");
    hFeSignal->SetYTitle("# counts");
        
        // Write convolution histograms in root files
        TFile* f = new TFile(fOutputName, "RECREATE");
        hFeElectrons->Write();
        hFeSignal->Write();
        f->Close();
    //}

    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    return 0;
}


