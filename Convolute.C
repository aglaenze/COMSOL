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

#include "Include/Utils.C"


/* Creates histograms that are the convolution of Fe55 spectrum histogram (number of primaries created by 1 photon) and secondaries histogram (number of secondaries created by one primary in the drift region, after amplification)
 */
using namespace std;

int Convolute(int modelNum, string gasName, vector<int> hvList) {
	
	
	const bool drawConvoluteSpectrum = false;
	
	
	time_t t0 = time(NULL);
	
	//gStyle->SetOptStat(0);
	
	const int electrodeNum = GetElectrodeNum(modelNum);
	if ((int)hvList.size() != electrodeNum-1) {cout << "Wrong hv input" << endl; return 0;}
	
	const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	Int_t num = GetNumberOfFiles(path, "signal");
	num = 1;
	
	Int_t nPrimaryTh = GetPrimary(gasName);
	
	TString fInputHV;
	
	for (int k = 0; k< (int)hvList.size(); k++) fInputHV += Form("-%d", hvList[k]);
	fInputHV += ".root";
	
	TString fSignalName = path+ "signal" + fInputHV;
	TString fOutputName = path+ "fe-spectrum-convoluted" + fInputHV;
	
	TFile* f = new TFile(fOutputName, "RECREATE");
	
	map <string, int, NoSorting> electrode;
	LoadElectrodeMap(modelNum, electrode);
	int padElectrode = 0;
	int driftElectrode = 0;
	
	map<string, int>::iterator it = electrode.begin();
	for (it=electrode.begin(); it!=electrode.end(); ++it) {
		//cout << it->first << " => " << it->second << '\n';
		if (it->first == "pad") padElectrode = it->second;
		else if (it->first == "drift") driftElectrode = it->second;
	}
	if (padElectrode == 0 || driftElectrode == 0) {cout << "Did not find drift or pad electrode" << endl; return 0;}
	
	// open input files
	TFile* fFe = new TFile(Form("rootFiles/%s/spectrum_Fe55.root", gasName.c_str()), "read");
	TH1F* hFe = (TH1F*)fFe->Get("hElectrons");
	Int_t nFe = hFe->GetEntries();
	
	TFile* fSignal = new TFile(fSignalName, "read");
	
	
	// Initialise the tree of avalanche size
	TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
	tAvalanche->SetBranchStatus("*",0);	// Select the branches to look at (do't include the vectors!!)
	tAvalanche->SetBranchStatus("amplificationElectrons",1);
	tAvalanche->SetBranchStatus("ionNum",1);
	tAvalanche->SetBranchStatus("ionBackNum",1);
	
	Int_t nAvalanche = tAvalanche->GetEntries();
	//Int_t nAmplification;
	vector <Int_t> *nWinnersVec = nullptr;
	vector <Int_t> gainVec = {};
	Int_t ni = 0, ionBackNum = 0;
	//tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
	tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVec);
	tAvalanche->SetBranchAddress("ionNum", &ni);
	tAvalanche->SetBranchAddress("ionBackNum", &ionBackNum);

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
			cout << "nAvalanche != nCharge" << endl;
			return 0;
		}
		cout << "Number of entries in tAvalanche = " << nAvalanche << endl;
		cout << "Number of entries in tCharge = " << nCharge << endl;
		cout << "Number of entries in hFe = " << nFe << endl;
		hFeCharge[k] = new TH1F(Form("hFeCharge_%d", k+2), "Number of induced charges with Fe source", nBins, 0, nBins );
		hFeCharge[k]->SetXTitle("# induced charges");
		hFeCharge[k]->SetYTitle("# counts");
	}
	
	//return 0;
	// Convolute the trees of the avalanche size and charges with Fe spectrum
	TH1F* hFeAmplification = new TH1F("hFeAmplification", "Number of avalanche electrons with Fe source", nBins, 0, nBins );
	TH1F* hFeIbf = new TH1F("hFeIbf", "IBF with Fe source", 10000, 0, 100 );
	TH1F* hFeIbfTotalCharge = new TH1F("hFeIbfTotalCharge", "Total induced charge IBF with Fe source", 10000, 0, 100 );
	TH1F* hFeIbfIonCharge = new TH1F("hFeIbfIonCharge", "Induced ion charge IBF with Fe source", 10000, 0, 100 );
	
	const int numberOfPhotons = 5000;
	for (unsigned int i = 0; i < numberOfPhotons; ++i) {
		if (i % (int(numberOfPhotons/10)) == 0) cout << i << "/" << numberOfPhotons << " photons" << endl;
		Int_t nPrim = hFe->GetRandom();
		//cout << "\nNprim = " << nPrim << endl;
		Double_t gtot = 0, itot = 0, ibntot = 0;
		Double_t ctot[electrodeNum], ctotIon[electrodeNum];
		for (int k = 0; k<electrodeNum; k++) {ctot[k] = 0; ctotIon[k] = 0;}
		for (Int_t j = 0; j < nPrim; ++j) {
			int r = rand() % nAvalanche;
			tAvalanche->GetEntry(r);
			//gtot += nAmplification;
			gainVec = *nWinnersVec;
			gtot += gainVec[0];
			itot += ni;
			ibntot += ionBackNum;
			for (int k = 0; k<electrodeNum; k++) {
				tCharge[k]->GetEntry(r);
				ctot[k] += abs(totalInducedCharge[k]);
				ctotIon[k] += abs(ionInducedCharge[k]);
			}
		}
		hFeAmplification->Fill(gtot/nPrimaryTh);
		if (itot > 10 ) hFeIbf->Fill((double)ibntot/itot * 100.);
		for (int k = 0; k<electrodeNum; k++) {hFeCharge[k]->Fill(ctot[k]/nPrimaryTh);}
		if (ctot[padElectrode-2]>10){
			hFeIbfTotalCharge->Fill(ctot[driftElectrode-2]/ctot[padElectrode-2]*100.);
			hFeIbfIonCharge->Fill(ctotIon[driftElectrode-2]/ctotIon[padElectrode-2]*100.);
		}
	}
	//return 0;
	
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
	hFeIbf->Write();
	for (int k = 0; k < electrodeNum; k++) {hFeCharge[k]->Write();}
	hFeIbfTotalCharge->Write();
	hFeIbfIonCharge->Write();
	f->Close();
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	
	return 0;
}


