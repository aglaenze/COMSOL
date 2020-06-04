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



int Transparency() {
	
	//______________________
	// variables
	std::string gasName = "Ne"; // Ar-iC4H10 or Ne or Ar-CO2
	const int modelNum = 15;
	//____________________
	
	time_t t0 = time(NULL);
	gStyle->SetTitleFontSize(.06);
	gStyle->SetTitleSize(.06);
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetLabelSize(.05, "XY");
	gStyle->SetMarkerSize(0.3);
	gStyle->SetTextSize(0.05);
	
	const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	TString fileName = path + "fe-spectrum-convoluted-360-405-690-830-910-1230.root";
	
	std::map <std::string, int> electrode;
	LoadElectrodeMap(modelNum, electrode);
	
	int driftElectrode = electrode["drift"];
	int gem2upElectrode = electrode["GEM2 up"];
	int gem2downElectrode = electrode["GEM2 down"];
	int gem1upElectrode = electrode["GEM1 up"];
	int gem1downElectrode = electrode["GEM1 down"];
	int meshElectrode = electrode["mesh"];
	int padElectrode = electrode["pad"];
	
	std::vector <Int_t> hv = {};
	TObjArray* matches = TPRegexp( "-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?\\.root$" ).MatchS( fileName );
	for (int k = 0; k<matches->GetLast(); k++) {
		hv.push_back( (static_cast<TObjString*>(matches->At(k+1)))->GetString().Atoi() );
		std::cout << hv[k] << std::endl;
	}
	
	TCanvas* cv = new TCanvas("cv","cv", 1000, 800);
	//cv->Divide(2);
	cv->Divide(2, 2);
	
	TString fSignalName = path + "signal";
	for (int k = 0; k<hv.size(); k++) {fSignalName += Form("-%d", hv[k]);}
	fSignalName += ".root";
	std::cout << fSignalName << std::endl;
	
	TFile* fSignal = TFile::Open(fSignalName, "READ");
	
	/*
	 // To draw the transparency I need to get:
	 1) full detector: if there are charges at all in the end (pads)
	 2) GEM2: charges readout on GEM1 up ?
	 3) GEM1: IF there were charges "sensed" by GEM1 up, were there charges "sensed" by the mesh? --> plot gain x transparency
	 4) mesh: IF there were charges "sensed" by GEM1 down, how many  were sensed  by the pads --> plot gain x transparency
	 
	 */
	
	TTree* tChargeDrift = (TTree*)fSignal->Get(Form("tInducedCharge_%d", driftElectrode));
	TTree* tChargeGem2up = (TTree*)fSignal->Get(Form("tInducedCharge_%d", gem2upElectrode));
	TTree* tChargeGem2down = (TTree*)fSignal->Get(Form("tInducedCharge_%d", gem2downElectrode));
	TTree* tChargeGem1up = (TTree*)fSignal->Get(Form("tInducedCharge_%d", gem1upElectrode));
	TTree* tChargeGem1down = (TTree*)fSignal->Get(Form("tInducedCharge_%d", gem1downElectrode));
	TTree* tChargeMesh = (TTree*)fSignal->Get(Form("tInducedCharge_%d", meshElectrode));
	TTree* tChargePad = (TTree*)fSignal->Get(Form("tInducedCharge_%d", padElectrode));
	
	Double_t chargeDrift, chargeGem2up, chargeGem2down, chargeGem1up, chargeGem1down, chargeMesh, chargePad;
	tChargeDrift->SetBranchAddress("electronInducedCharge", &chargeDrift);
	tChargeGem2up->SetBranchAddress("electronInducedCharge", &chargeGem2up);
	tChargeGem2down->SetBranchAddress("electronInducedCharge", &chargeGem2down);
	tChargeGem1up->SetBranchAddress("electronInducedCharge", &chargeGem1up);
	tChargeGem1down->SetBranchAddress("electronInducedCharge", &chargeGem1down);
	tChargeMesh->SetBranchAddress("electronInducedCharge", &chargeMesh);
	tChargePad->SetBranchAddress("electronInducedCharge", &chargePad);
	
	
	int nChargePad = tChargePad->GetEntries();
	int nChargeDrift = tChargeDrift->GetEntries();
	if (nChargePad != nChargeDrift) {
		std::cout << "nChargePad != nChargeDrift" << std::endl;
		return 0;
	}
	int nCharge = nChargeDrift;
	
	int winners = 0, winners2 = 0, winners3 = 0, winners4 = 0;
	int ntotal3 = 0, ntotal4 = 0;
	TH1F* hTransparency = new TH1F("hTransparency", "hTransparency", 4, -1, 3);
	TH1F* hTransparencyGem2 = new TH1F("hTransparencyGem2", "hTransparencyGem2", 4, -1, 3);
	TH1F* hTransparencyGem1 = new TH1F("hTransparencyGem1", "hTransparencyGem1", 4, -1, 3);
	TH1F* hTransparencyMesh = new TH1F("hTransparencyMesh", "hTransparencyMesh", 4, -1, 3);
	//TH1F* hTransparencyGem1 = new TH1F("hTransparencyGem1", "hTransparencyGem1", 100, 0, 6);
	//TH1F* hTransparencyMesh = new TH1F("hTransparencyMesh", "hTransparencyMesh", 100, 0, 100);
	for (int k = 0; k<nCharge; k++) {
		tChargeGem2up->GetEntry(k);
		tChargeGem2down->GetEntry(k);
		tChargeGem1up->GetEntry(k);
		tChargeGem1down->GetEntry(k);
		tChargeMesh->GetEntry(k);
		tChargePad->GetEntry(k);
		// 1) full detector: if there are charges at all in the end (pads)
		if (abs(chargeDrift)>0) {
			if (abs(chargePad)>1) { hTransparency->Fill(1); winners++;}
			else {hTransparency->Fill(0);}
		}
		// 2) GEM2: charges readout on GEM1 up ?
		if (abs(chargeDrift)>0) {
			//std::cout << "OK " << std::endl;
			if (abs(chargeGem1up)>1) {hTransparencyGem2->Fill(1); winners2++;}
			else {hTransparencyGem2->Fill(0);}
		}
		//else std::cout << "abs(chargeGem2up) = " << abs(chargeGem2up) << std::endl;
		// 3) GEM1: IF there were charges "sensed" by GEM1 up, were there charges "sensed" by the mesh? --> plot gain x transparency
		if (abs(chargeGem1up)>1) {
			ntotal3++;
			//std::cout << "abs(chargeGem1up) = " << abs(chargeGem1up) << std::endl;
			if (abs(chargeMesh)>1) {hTransparencyGem1->Fill(1); winners3++;}
			else hTransparencyGem1->Fill(0);
			//hTransparencyGem1->Fill(abs(chargeGem1down/chargeGem1up));
		}
		//4) mesh: IF there were charges "sensed" by GEM1 down, how many  were sensed  by the pads --> plot gain x transparency
		if (abs(chargeGem1down)>1) {
			ntotal4++;
			if (abs(chargePad)>1) {hTransparencyMesh->Fill(1); winners4++;}
			else hTransparencyMesh->Fill(0);
			//hTransparencyMesh->Fill(abs(chargePad/chargeGem1down));
		}
	}
	
	// Draw the full transparency of the detector
	cv->cd(1);
	hTransparency->SetTitle("Full transparency");
	hTransparency->Scale(1./nCharge);
	hTransparency->SetMaximum(1.2);
	hTransparency->Draw("hist");
	// Write text with the value of transparency
	Double_t transp = (double)winners/nCharge;
	
	TText* txttr = new TText(0.2,0.9*hTransparency->GetMaximum(),Form("Transparency = %.1f %s", transp*100, "%"));
	txttr->Draw();
	
	// Draw the transparency of GEM2
	cv->cd(2);
	hTransparencyGem2->SetTitle("GEM2 transparency");
	hTransparencyGem2->Scale(1./nCharge);
	hTransparencyGem2->SetMaximum(1.2);
	hTransparencyGem2->Draw("hist");
	// Write text with the value of transparency
	Double_t transp2 = (double)winners2/nCharge;
	TText* txttr2 = new TText(0.2,0.9*hTransparencyGem2->GetMaximum(),Form("Transparency = %.1f %s", transp2*100, "%"));
	txttr2->Draw();
	
	// Draw the transparency of GEM1
	cv->cd(3);
	hTransparencyGem1->SetTitle("GEM1 transparency");
	hTransparencyGem1->Scale(1./nCharge);
	hTransparencyGem1->SetMaximum(1.2);
	hTransparencyGem1->Draw("hist");
	
	// Write text with the value of transparency
	Double_t transp3 = (double)winners3/ntotal3;
	TText* txttr3 = new TText(0.2,0.9*hTransparencyGem1->GetMaximum(),Form("Transparency = %.1f %s", transp3*100, "%"));
	txttr3->Draw();
	
	
	// Draw the transparency of GEM2
	cv->cd(4);
	hTransparencyMesh->SetTitle("Mesh transparency");
	hTransparencyMesh->Scale(1./nCharge);
	hTransparencyMesh->SetMaximum(1.2);
	hTransparencyMesh->Draw("hist");
	// Write text with the value of transparency
	Double_t transp4 = (double)winners4/ntotal4;
	TText* txttr4 = new TText(0.2,0.9*hTransparencyMesh->GetMaximum(),Form("Transparency = %.1f %s", transp4*100, "%"));
	txttr4->Draw();
	
	
	
	TString outputName = Form("Figures/model%d/transparency-%s", modelNum, gasName.c_str());
	for (int k = 0; k<hv.size(); k++) {outputName += Form("-%d", hv[k]);}
	outputName+=".pdf";
	cv->SaveAs(outputName);
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	
	return 0;
}


