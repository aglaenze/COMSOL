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

int AddSignalTrees(int modelNum, std::string gasName, std::vector<int> hvList, bool computeIbf, bool feSource) {
	
	
	time_t t0 = time(NULL);
	
	const int electrodeNum = GetElectrodeNum(modelNum);
	if ((int)hvList.size() != electrodeNum-1) {std::cout << "Wrong HV input" << std::endl; return 0;}
	
	TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	
	TString filename = "signal";
	if (feSource) filename = "fe-signal";
	if (!computeIbf) filename += "-noibf";
	for (int k = 0; k< (int)hvList.size(); k++) filename += Form("-%d", hvList[k]);
	std::cout << "processing " << filename << "* (model " << modelNum << ")" << std::endl;
	
	int numberOfFiles = GetNumberOfFiles(path, filename+"-");
	std::cout << "number of files = " << numberOfFiles << std::endl;
	if (numberOfFiles == 0) {std::cout << "no files to add, terminating" << std::endl; return 0;}
	
	TString outputName = path + filename+".root";
	TFile* fOut = new TFile(outputName, "RECREATE");
	
	// Initialisation of the new TTree tAvalanche
	TTree *tAvalanche = new TTree("tAvalanche","Gain");
	Int_t nWinners = 0, ne2 = 0;
	Int_t ni = 0, ionBackNum = 0;
	std::vector<float> electronStartPoints = {}, electronEndPoints = {};
	std::vector<float> ionStartPoints = {}, ionEndPoints = {}, ionEndPointsX = {}, ionEndPointsY = {};
	tAvalanche->Branch("amplificationElectrons", &nWinners, "amplificationElectrons/I");
	tAvalanche->Branch("avalancheSize", &ne2, "avalancheSize/I");
	if (computeIbf) {
		tAvalanche->Branch("ionNum", &ni, "ibfRatio/I");
		tAvalanche->Branch("ionBackNum", &ionBackNum, "ionBackNum/I");
		tAvalanche->Branch("electronStartPoints", &electronStartPoints);
		tAvalanche->Branch("electronEndPoints", &electronEndPoints);
		tAvalanche->Branch("ionStartPoints", &ionStartPoints);
		tAvalanche->Branch("ionEndPoints", &ionEndPoints);
		tAvalanche->Branch("ionEndPointsX", &ionEndPointsX);
		tAvalanche->Branch("ionEndPointsY", &ionEndPointsY);
	}
	
	std::vector<float> *electronStartPointsInput = 0, *electronEndPointsInput = 0;
	std::vector<float> *ionStartPointsInput = 0, *ionEndPointsInput = 0, *ionEndPointsInputX = 0, *ionEndPointsInputY = 0;
	
	// On commence par remplir le TTree tAvalanche
	for (int i = 1; i<numberOfFiles+1; i++) {
		//for (int i = 2; i<3; i++) {
		TString inputName = path + filename + Form("-%d.root", i);
		std::cout << "including " << inputName << std::endl;
		TFile* fIn = TFile::Open(inputName, "READ");
		TTree* tAvalancheIn = (TTree*)fIn->Get("tAvalanche");
		tAvalancheIn->SetBranchAddress("amplificationElectrons", &nWinners);
		tAvalancheIn->SetBranchAddress("avalancheSize", &ne2);
		if (computeIbf) {
			tAvalancheIn->SetBranchAddress("ionNum", &ni);
			tAvalancheIn->SetBranchAddress("ionBackNum", &ionBackNum);
			
			tAvalancheIn->SetBranchAddress("electronStartPoints", &electronStartPointsInput);
			tAvalancheIn->SetBranchAddress("electronEndPoints", &electronEndPointsInput);
			tAvalancheIn->SetBranchAddress("ionStartPoints", &ionStartPointsInput);
			tAvalancheIn->SetBranchAddress("ionEndPoints", &ionEndPointsInput);
			tAvalancheIn->SetBranchAddress("ionEndPointsX", &ionEndPointsInputX);
			tAvalancheIn->SetBranchAddress("ionEndPointsY", &ionEndPointsInputY);
		}
		
		int nIn = tAvalancheIn->GetEntries();
		
		for (int l = 0; l<nIn; l++) {
			tAvalancheIn->GetEntry(l);
			if (!computeIbf) {tAvalanche->Fill(); continue;}
			
			electronStartPoints = *electronStartPointsInput;
			electronEndPoints = *electronEndPointsInput;
			ionStartPoints = *ionStartPointsInput;
			ionEndPoints = *ionEndPointsInput;
			ionEndPointsX = *ionEndPointsInputX;
			ionEndPointsY = *ionEndPointsInputY;
			
			tAvalanche->Fill();
			electronStartPoints.clear();
			electronEndPoints.clear();
			ionStartPoints.clear();
			ionEndPoints.clear();
			ionEndPointsX.clear();
			ionEndPointsY.clear();
		}
		
	}
	fOut->cd();
	tAvalanche->Write();
	
	if (!computeIbf) {
		fOut->Close();
		time_t t1 = time(NULL);
		PrintTime(t0, t1);
		return 0;
	}
	
	//return 0;
	
	for (int k = 0; k < electrodeNum; k++) {
		std::cout << "looping over electrode " << k+2 << "..." << std::endl;
		// Initialisation des TTrees tCharge
		
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
			
			TTree* tChargeIn = (TTree*)fIn->Get(Form("tInducedCharge_%d",k+2));
			tChargeIn->SetBranchAddress("totalInducedCharge", &tic);
			tChargeIn->SetBranchAddress("electronInducedCharge", &eic);
			tChargeIn->SetBranchAddress("ionInducedCharge", &iic);
			
			
			int nChargeIn = tChargeIn->GetEntries();
			for (int l = 0; l<nChargeIn; l++) {
				tChargeIn->GetEntry(l);
				tCharge->Fill();
			}
			fIn->Close();
		}
		fOut->cd();
		tCharge->Write("", TObject::kOverwrite);
		
	}
	fOut->Close();
	
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	return 0;
	
}


