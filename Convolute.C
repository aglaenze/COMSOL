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
	
	const int electrodeNum = GetElectrodeNum(modelNum);
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
	
	std::map <std::string, int> electrode;
	LoadElectrodeMap(modelNum, electrode);
	int padElectrode = electrode["pad"];
	int driftElectrode = electrode["drift"];
	
	// open input files
	TFile* fFe = new TFile(Form("rootFiles/%s/spectrum_Fe55.root", gasName.c_str()), "read");
	TH1F* hFe = (TH1F*)fFe->Get("hElectrons");
	Int_t nFe = hFe->GetEntries();
	
	TFile* fSignal = new TFile(fSignalName, "read");
	
	// Initialise the tree of avalanche size
	TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
	Int_t nAvalanche = tAvalanche->GetEntries();
	Int_t nAmplification;
	tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
	const Int_t nBins = int(tAvalanche->GetMaximum("amplificationElectrons"));
	
	// Initialise the trees of induced charges
	TH1F* hFeCharge[electrodeNum];
	TTree* tCharge[electrodeNum];
	Double_t totalInducedCharge[electrodeNum], ionInducedCharge[electrodeNum];
	for (int k = 0; k<electrodeNum; k++) {
		tCharge[k] = (TTree*) fSignal->Get(Form("tInducedCharge_%d", k+2));
		tCharge[k]->SetBranchAddress("totalInducedCharge", &totalInducedCharge[k]);
		tCharge[k]->SetBranchAddress("ionInducedCharge", &ionInducedCharge[k]);
		Int_t nCharge = tCharge[k]->GetEntries();
		if (nAvalanche != nCharge) {
			std::cout << "nAvalanche != nCharge" << std::endl;
			return 0;
		}
		std::cout << "Number of entries in tAvalanche = " << nAvalanche << std::endl;
		std::cout << "Number of entries in tCharge = " << nCharge << std::endl;
		std::cout << "Number of entries in hFe = " << nFe << std::endl;
		hFeCharge[k] = new TH1F(Form("hFeCharge_%d", k+2), "Number of induced charges with Fe source", nBins, 0, nBins );
		hFeCharge[k]->SetXTitle("# induced charges");
		hFeCharge[k]->SetYTitle("# counts");
	}
	
	/*
	 // New temporary TTree
	 TTree *tChargeConvoluted = new TTree("tChargeConvoluted","Convoluted Charges");
	 Double_t padIonCharge = 0, padTotalCharge = 0, driftIonCharge = 0, driftTotalCharge = 0;
	 tChargeConvoluted->Branch("padIonCharge", &padIonCharge, "padIonCharge/D");
	 tChargeConvoluted->Branch("padTotalCharge", &padTotalCharge, "padTotalCharge/D");
	 tChargeConvoluted->Branch("driftIonCharge", &driftIonCharge, "driftIonCharge/D");
	 tChargeConvoluted->Branch("driftTotalCharge", &driftTotalCharge, "driftTotalCharge/D");
	 */
	
	
	// Convolute the trees of the avalanche size and charges with Fe spectrum
	TH1F* hFeAmplification = new TH1F("hFeAmplification", "Number of avalanche electrons with Fe source", nBins, 0, nBins );
	TH1F* hFeIbfTotal = new TH1F("hFeIbfTotal", "Total induced charge IBF with Fe source", 10000, 0, 100 );
	TH1F* hFeIbfIon = new TH1F("hFeIbfIon", "Induced ion charge IBF with Fe source", 10000, 0, 100 );
	for (unsigned int i = 0; i < 10000; ++i) {
		Int_t nPrim = hFe->GetRandom();
		//std::cout << "\nNprim = " << nPrim << std::endl;
		Double_t gtot = 0;
		Double_t ctot[electrodeNum], ctotIon[electrodeNum];
		for (int k = 0; k<electrodeNum; k++) {ctot[k] = 0; ctotIon[k] = 0;}
		for (unsigned int j = 0; j < nPrim; ++j) {
			int r = rand() % nAvalanche;
			tAvalanche->GetEntry(r);
			gtot += nAmplification;
			for (int k = 0; k<electrodeNum; k++) {
				tCharge[k]->GetEntry(r);
				ctot[k] += abs(totalInducedCharge[k]);
				ctotIon[k] += abs(ionInducedCharge[k]);
			}
		}
		hFeAmplification->Fill(gtot/nPrimaryTh);
		for (int k = 0; k<electrodeNum; k++) {hFeCharge[k]->Fill(ctot[k]/nPrimaryTh);}
		if (ctot[padElectrode-2]>10){
			hFeIbfTotal->Fill(ctot[driftElectrode-2]/ctot[padElectrode-2]*100.);
			hFeIbfIon->Fill(ctotIon[driftElectrode-2]/ctotIon[padElectrode-2]*100.);
		}
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
	for (int k = 0; k<electrodeNum; k++) {hFeCharge[k]->Write();}
	hFeIbfTotal->Write();
	hFeIbfIon->Write();
	f->Close();
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	
	return 0;
}


