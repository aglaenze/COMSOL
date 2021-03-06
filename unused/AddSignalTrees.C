#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "_Utils.C"

/*
 I simulate root files with 100 events, and if I want to have them as one file, I have to add current
 */

// If I also want to keep signal trees (unused now)

int AddSignalTrees(int modelNum, std::string gasName, std::vector<int> hvList) {
	
	
	time_t t0 = time(NULL);
	
	const int electrodeNum = GetElectrodeNum(modelNum);
	if (hvList.size() != electrodeNum-1) {std::cout << "Wrong hv input" << std::endl; return 0;}
	
	TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	
	TString filename = "signal";
	for (int k = 0; k< hvList.size(); k++) filename += Form("-%d", hvList[k]);
	std::cout << "processing " << filename << " (model " << modelNum << ")" << std::endl;
	//return 0;
	
	
	int numberOfFiles = GetNumberOfFiles(path, filename+"-");
	std::cout << "number of files = " << numberOfFiles << std::endl;
	TString outputName = path + filename+".root";
	//return 0;
	TFile* fOut = new TFile(outputName, "RECREATE");
	
	
	// Initialisation of the new TTree tAvalanche
	TTree *tAvalanche = new TTree("tAvalanche","Gain");
	Int_t nWinners = 0, ne2 = 0;
	//Double_t ibfRatio = 0.;
	int ni = 0, ionBackNum = 0;
	tAvalanche->Branch("amplificationElectrons", &nWinners, "amplificationElectrons/I");
	tAvalanche->Branch("avalancheSize", &ne2, "avalancheSize/I");
	//tAvalanche->Branch("ibfRatio", &ibfRatio, "ibfRatio/D");
	tAvalanche->Branch("ionNum", &ni, "ibfRatio/I");
	tAvalanche->Branch("ionBackNum", &ionBackNum, "ionBackNum/I");
	
	// On commence par remplir le TTree tAvalanche
	for (int i = 1; i<numberOfFiles+1; i++) {
		//for (int i = 2; i<3; i++) {
		TString inputName = path + filename + Form("-%d.root", i);
		TFile* fIn = TFile::Open(inputName, "READ");
		TTree* tAvalancheIn = (TTree*)fIn->Get("tAvalanche");
		tAvalancheIn->SetBranchAddress("amplificationElectrons", &nWinners);
		tAvalancheIn->SetBranchAddress("avalancheSize", &ne2);
		//tAvalancheIn->SetBranchAddress("ibfRatio", &ibfRatio);
		tAvalancheIn->SetBranchAddress("ionNum", &ni);
		tAvalancheIn->SetBranchAddress("ionBackNum", &ionBackNum);
		
		int nIn = tAvalancheIn->GetEntries();
		
		for (int l = 0; l<nIn; l++) {
			tAvalancheIn->GetEntry(l);
			//if (ibfRatio > 1) {std::cout << "problem!" << std::endl; return 0;}
			tAvalanche->Fill();
		}
		
	}
	fOut->cd();
	tAvalanche->Write();
	
	//return 0;
	
	for (int k = 0; k < electrodeNum; k++) {
		std::cout << "looping over electrode " << k+2 << "..." << std::endl;
		// Initialisation des TTrees tSignal et tCharge
		TTree *tSignal = new TTree(Form("tSignal_%d",k+2),"Currents");
		Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
		Double_t fctInt = 0, fceInt = 0, fciInt = 0;
		tSignal->Branch("time", &ft, "time/D");
		tSignal->Branch("totalCurrent", &fct, "totalCurrent/D");
		tSignal->Branch("electronCurrent", &fce, "electronCurrent/D");
		tSignal->Branch("ionCurrent", &fci, "ionCurrent/D");
		tSignal->Branch("totalCurrentInt",&fctInt,"totalCurrentInt/D");     // Integral
		tSignal->Branch("electronCurrentInt",&fceInt,"electronCurrentInt/D");
		tSignal->Branch("ionCurrentInt",&fciInt,"ionCurrentInt/D");
		
		TTree *tCharge = new TTree(Form("tInducedCharge_%d",k+2),"InducedCharge");
		Double_t tic = 0., eic = 0., iic = 0.;
		tCharge->Branch("totalInducedCharge", &tic, "totalInducedCharge/D");
		tCharge->Branch("electronInducedCharge", &eic, "electronInducedCharge/D");
		tCharge->Branch("ionInducedCharge", &iic, "ionInducedCharge/D");
		
		// Maintenant on va remplir les TTrees tSignal et tCharge
		Double_t tFin = 0;
		Double_t timeStep = 0;
		for (int i = 1; i<numberOfFiles+1; i++) {
			TString inputName = path + filename + Form("-%d.root", i);
			TFile* fIn = TFile::Open(inputName, "READ");
			
			TTree* tSignalIn = (TTree*)fIn->Get(Form("tSignal_%d",k+2));
			Double_t timeIn;
			tSignalIn->SetBranchAddress("time", &timeIn);
			tSignalIn->SetBranchAddress("totalCurrent", &fct);
			tSignalIn->SetBranchAddress("electronCurrent", &fce);
			tSignalIn->SetBranchAddress("ionCurrent", &fci);
			tSignalIn->SetBranchAddress("totalCurrentInt", &fctInt);
			tSignalIn->SetBranchAddress("electronCurrentInt", &fceInt);
			tSignalIn->SetBranchAddress("ionCurrentInt", &fciInt);
			
			TTree* tChargeIn = (TTree*)fIn->Get(Form("tInducedCharge_%d",k+2));
			tChargeIn->SetBranchAddress("totalInducedCharge", &tic);
			tChargeIn->SetBranchAddress("electronInducedCharge", &eic);
			tChargeIn->SetBranchAddress("ionInducedCharge", &iic);
			
			// obtention de la time step (si on ne regarde pas tSignal on s'en fiche)
			tSignalIn->GetEntry(3); // nombre arbitraire
			Double_t t0 = timeIn;
			tSignalIn->GetEntry(4);
			Double_t t1 = timeIn;
			timeStep = t1-t0;
			//std::cout << "timeStep = " << timeStep << std::endl;
			
			int nSignalIn = tSignalIn->GetEntries();
			for (int l = 0; l<nSignalIn; l++) {
				tSignalIn->GetEntry(l);
				ft = tFin+timeIn;
				tSignal->Fill();
			}
			tFin = ft + timeStep;
			
			int nChargeIn = tChargeIn->GetEntries();
			for (int l = 0; l<nChargeIn; l++) {
				tChargeIn->GetEntry(l);
				tCharge->Fill();
			}
			fIn->Close();
		}
		fOut->cd();
		tSignal->Write("", TObject::kOverwrite);
		tCharge->Write("", TObject::kOverwrite);
		
		std::cout << "tMax = " << tSignal->GetMaximum("time") << std::endl;

		
	}
	fOut->Close();
	
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	return 0;
	
}


