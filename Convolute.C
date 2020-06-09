#include <iostream>
#include <fstream>
#include <cmath>
#include <dirent.h>

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "_Utils.C"

#include <cmath>

/* Creates histograms that are the convolution of Fe55 spectrum histogram (number of primaries created by 1 photon) and secondaries histogram (number of secondaries created by one primary in the drift region, after amplification)
 */


int Convolute(int modelNum, std::string gasName, std::vector<int> hvList) {
    

    const bool drawConvoluteSpectrum = false;
    
    
    time_t t0 = time(NULL);
    
    gStyle->SetOptStat(0);

	int electrodeNum = GetElectrodeNum(modelNum);
	if (hvList.size() != electrodeNum-1) {std::cout << "Wrong hv input" << std::endl; return 0;}
     
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    Int_t num = GetNumberOfFiles(path, "signal");
    num = 1;
    
    Int_t nPrimaryTh = GetPrimary(gasName);
    
    TString fInputHV;
	
	for (int k = 0; k< hvList.size(); k++) fInputHV += Form("-%d", hvList[k]);
	fInputHV += ".root";
    
	TString fSignalName = path+ "signal" + fInputHV;
	TString fOutputName = path+ "fe-spectrum-convoluted" + fInputHV;

    TFile* f = new TFile(fOutputName, "RECREATE");
    
    //int readoutElectrode = electrode["mesh"];   // if readout electrode = pad, do not foget the - sign
    
    // open input files
    TFile* fFe = new TFile(Form("rootFiles/%s/spectrum_Fe55.root", gasName.c_str()), "read");
    TH1F* hFe = (TH1F*)fFe->Get("hElectrons");
    Int_t nFe = hFe->GetEntries();
    
    TFile* fSignal = new TFile(fSignalName, "read");
    TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
    Int_t nAvalanche = tAvalanche->GetEntries();
    
    Int_t nAmplification;
    tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
    const Int_t nBins = int(tAvalanche->GetMaximum("amplificationElectrons"));
	
	// New TTree

	//tChargeConvoluted->Branch("test", ionCharge, "test/D");
	//return 0;
	
	TTree *tChargeConvoluted = new TTree("tChargeConvoluted","Convoluted Charges");
    
    // First convolute the trees of induced charges with Fe spectrum
    for (int k = 0; k<electrodeNum; k++) {
        TTree* tCharge = (TTree*) fSignal->Get(Form("tInducedCharge_%d", k+2));
        Double_t totalInducedCharge;
        tCharge->SetBranchAddress("totalInducedCharge", &totalInducedCharge);
        Int_t nCharge = tCharge->GetEntries();
        if (nAvalanche != nCharge) {
            std::cout << "nAvalanche != nCharge" << std::endl;
            return 0;
        }
        std::cout << "Number of entries in tAvalanche = " << nAvalanche << std::endl;
        std::cout << "Number of entries in tCharge = " << nCharge << std::endl;
        std::cout << "Number of entries in hFe = " << nFe << std::endl;
        TH1F* hFeCharge = new TH1F(Form("hFeCharge_%d", k+2), "Number of induced charges with Fe source", nBins, 0, nBins );
        
		Double_t totalCharge;
		TBranch* b = tChargeConvoluted->Branch(Form("totalInducedCharge_%d", k+2), &totalCharge, Form("totalInducedCharge_%d/D", k+2));
		
        for (unsigned int i = 0; i < 10000; ++i) {
            Int_t nPrim = hFe->GetRandom();
            //std::cout << "\nNprim = " << nPrim << std::endl;
            Double_t gtot = 0;
            for (unsigned int j = 0; j < nPrim; ++j) {
                int r = rand() % nCharge;
                tCharge->GetEntry(r);
				// the following line is a test
				//if (abs(totalInducedCharge) == 0) {j--; continue;}
                gtot += abs(totalInducedCharge);
            }
            hFeCharge->Fill(gtot/nPrimaryTh);
			totalCharge = gtot/nPrimaryTh;
			//b->Fill();
			
        }
        hFeCharge->SetXTitle("# induced charges");
        hFeCharge->SetYTitle("# counts");
        // Write convolution histograms in root files
        f->cd();
        hFeCharge->Write();
    }
	//tChargeConvoluted->Write();

    //return 0;
    // Second convolute the trees of the avalanche size with Fe spectrum
    TH1F* hFeAmplification = new TH1F("hFeAmplification", "Number of avalanche electrons with Fe source", nBins, 0, nBins );
    for (unsigned int i = 0; i < 10000; ++i) {
        Int_t nPrim = hFe->GetRandom();
        //std::cout << "\nNprim = " << nPrim << std::endl;
        Double_t gtot = 0;
        for (unsigned int j = 0; j < nPrim; ++j) {
            int r = rand() % nAvalanche;
            tAvalanche->GetEntry(r);
			// the following line is a test
			if (nAmplification < 1) {j--; continue;}
            gtot += nAmplification;
        }
        hFeAmplification->Fill(gtot/nPrimaryTh);
    }
        
    if (drawConvoluteSpectrum) {
        TCanvas* c1 = new TCanvas();
        c1->cd();
            
        hFeAmplification->SetXTitle("# amplification electrons");
        hFeAmplification->SetYTitle("# counts");
        hFeAmplification->Draw("hist");
        
        c1->SaveAs(Form("Figures/model%d/Convolution.pdf", modelNum));
    }
        
    // Write convolution histogram in root files
    f->cd();
    hFeAmplification->Write();
    f->Close();
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    return 0;
}


