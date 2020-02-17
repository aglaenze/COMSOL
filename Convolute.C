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
    const int modelNum = 9;
    //____________________
    
    
    time_t t0 = time(NULL);
    
    gStyle->SetOptStat(0);

    TFile fFe(Form("rootFiles/%s/spectrum_Fe55.root", gasName.c_str()));
    TH1F* hFe = (TH1F*)fFe.Get("hElectrons");
     
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    Int_t num = GetNumberOfFiles(path, "gain");
    
    Int_t nPrimaryTh = GetPrimary(gasName);
    //Int_t nPrimaryTh = 600;
    
    TString fGainName, fOutputName;
    
    for (unsigned int k = 0; k < num; ++k) {
        Int_t hvMm = 0, hvDmDown = 0, hvDmUp = 0, hvDrift = 0;
        if (modelNum == 1) {
            hvMm = 340+20*k;
            hvDrift = 540+20*k;
            fGainName = path+Form("gain-%d-%d", hvMm, hvDrift);
            fOutputName = Form("Fe-spectrum-convoluted-%d-%d.root", hvMm, hvDrift);
        }
        else if (modelNum > 6 && modelNum < 10) {
            hvDmDown = 300;
            hvDmUp = 600;
            hvDrift = 800;
            fGainName = path+Form("gain-%d-%d-%d", hvDmDown, hvDmUp, hvDrift);
            fOutputName = Form("Fe-spectrum-convoluted-%d-%d-%d.root", hvDmDown, hvDmUp, hvDrift);
        }
        else if (modelNum == 10) {
            hvMm = 300;
            hvDmDown = 300;
            hvDmUp = 600;
            hvDrift = 800;
            fGainName = path+Form("gain-%d-%d-%d", hvDmDown, hvDmUp, hvDrift);
            fOutputName = Form("Fe-spectrum-convoluted-%d-%d-%d.root", hvDmDown, hvDmUp, hvDrift);
        }

        TFile fGain(fGainName + ".root");
        TH1F* hGain = (TH1F*) fGain.Get("hElectrons");
        //hGain->Rebin(15);
        
        Int_t nGain = hGain->GetEntries();
        Int_t nFe = hFe->GetEntries();
        std::cout << "Number of entries in hGain = " << nGain << std::endl;
        std::cout << "Number of entries in hFe = " << nFe << std::endl;

        const Int_t nBins = int(nGain/4);
        Int_t iBinMax = hGain->GetMaximumBin(); // Return location of bin with maximum value in the range
        Double_t xMax = hGain->GetXaxis()->GetBinCenter(iBinMax);
        //TH1F* hFeElectrons = new TH1F("hFeElectrons", "Number of secondary electrons with Fe source", nBins, 0, 0.02*TMath::Exp(0.0352*hvMm) );
        TH1F* hFeElectrons = new TH1F("hFeElectrons", "Number of secondary electrons with Fe source", nBins, 0, nBins*20 );
        
        std::cout << "maximum bin = " << hGain->GetMaximumBin() << std::endl;
        
        for (unsigned int i = 0; i < 100000; ++i) {
            Int_t nPrim = hFe->GetRandom();
            //std::cout << "\nNprim = " << nPrim << std::endl;
            Double_t gtot = 0;
            for (unsigned int j = 0; j < nPrim; ++j) {
                Double_t gain = hGain->GetRandom();
                //std::cout << "gainBin = " << gainBin << std::endl;
                //std::cout << "gain = " << gain << std::endl;
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
        TFile* f = new TFile(path+fOutputName, "RECREATE");
        hFeElectrons->Write();
        f->Close();
    }

    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    return 0;
}


