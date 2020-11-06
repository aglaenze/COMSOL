#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>

#include <TCanvas.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>

#include "Include/Functions.C"
#include "Include/Transparency.C"



int Spectrum(int modelNum, std::string gasName, std::vector<int> hvList) {
	
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
	
	int electrodeNum = GetElectrodeNum(modelNum);
	if ((int)hvList.size() != electrodeNum-1) {std::cout << "Wrong hv input" << std::endl; return 0;}
	
	//Load geometry parameters
	double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
	
	// input and output files
	const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	TString fileName = path + "fe-spectrum-convoluted";
	TString fSignalName = path + "fe-signal-noibf";
	TString outputName = Form("Figures/model%d/spectra-%s", modelNum, gasName.c_str());
	for (int k = 0; k< (int)hvList.size(); k++) {
		fileName += Form("-%d", hvList[k]);
		fSignalName += Form("-%d", hvList[k]);
		outputName += Form("-%d", hvList[k]);
	}
	fileName += ".root";
	fSignalName += ".root";
	outputName+=".pdf";
	
	std::cout << fSignalName << std::endl;
	
	TCanvas* cv = new TCanvas("cv","cv", 1200, 1000);
	//cv->Divide(2);
	
	// First draw convoluted spectrum
	TFile* fConvoluted = TFile::Open(fileName, "READ");
	TH1F* hFeAmplification = (TH1F*)fConvoluted->Get("hFeAmplification");
	while (hFeAmplification->GetMaximum() < 100) hFeAmplification->Rebin(2);
	
	hFeAmplification->Scale(1/hFeAmplification->GetMaximum());
	hFeAmplification->SetMaximum(1.2);
	//hFeAmplification->GetXaxis()->SetRangeUser(2, 10000);
	hFeAmplification->SetLineColor(kBlue);
	hFeAmplification->SetTitle("Gain convoluted with Fe source");
	
	TF1* f = GetFitCurve(hFeAmplification);
	f->SetLineColor(kBlue);
	
	Int_t iBinMax = hFeAmplification->GetMaximumBin();
	Double_t xMax = hFeAmplification->GetXaxis()->GetBinCenter( iBinMax );
	hFeAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
	cv->cd(1);
	//hFeAmplification->Draw("hist");
	//f->Draw("same");
	
	
	// Write gain, sigma and res = sigma/gain on the plot
	double xPos = f->GetParameter(0) + 2.5*f->GetParameter(1);
	TLatex* txtGain = new TLatex(xPos, 0.9, Form("Gain = %.1f #pm %.1f ", f->GetParameter(0), f->GetParError(0)));
	TLatex* txtSigma = new TLatex(xPos, 0.8, Form("Sigma = %.1f #pm %.1f ", f->GetParameter(1), f->GetParError(1)));
	Double_t resolution = f->GetParameter(1)/f->GetParameter(0);
	Double_t resolutionError = resolution * TMath::Sqrt( Square(f->GetParError(0)/f->GetParameter(0)) + Square(f->GetParError(1)/f->GetParameter(1)) );
	std::string percent = "%";
	TLatex* txtRes = new TLatex(xPos, 0.7, Form("Sigma/Mean = %.1f #pm %.1f %s", resolution*100, resolutionError*100, percent.c_str()));
	
	//txtGain->Draw("same"); txtSigma->Draw("same"); txtRes->Draw("same");
	
	// Second draw the gain with Fe source simulated
	TFile* fSignal = TFile::Open(fSignalName, "READ");
	TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
	Int_t nAvalanche = tAvalanche->GetEntries();
	/*
	 Int_t nAmplification;
	 tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
	 */
	vector<int> *neAvalVecIn = 0, *nWinnersVecIn = 0;
	tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVecIn);
	tAvalanche->SetBranchAddress("avalancheSize", &neAvalVecIn);
	vector<int> neAvalVec = {}, nWinnersVec = {};
	
	Int_t nPrimaryTh = GetPrimary(gasName);
	
	// Create histogram of amplification electrons
	//Int_t elMax = tAvalanche->GetMaximum("amplificationElectrons");
	//Int_t elMax = int(xMax*3);
	Int_t elMax = 200000;
	cout << "\n\n el max " << elMax << endl;
	TH1F* hAmplification = new TH1F("hAmplification", "Number of amplification electrons", int(elMax/1000), 0, elMax);
	for (int k = 0; k< (int)tAvalanche->GetEntries(); k++) {
		tAvalanche->GetEntry(k);
		nWinnersVec = *nWinnersVecIn;
		neAvalVec = *neAvalVecIn;
		int nAmplification = 0;
		//cout << nWinnersVec.size() << endl;
		for (int l = 0; l< (int)nWinnersVec.size(); l++) {nAmplification += nWinnersVec[l];}
		cout << nAmplification << endl;
		double normAmplicationEl = (double)nAmplification/(double)nPrimaryTh;
		if (nAmplification > nPrimaryTh) hAmplification->Fill(nAmplification);
		//cout << (double)nAmplification/nPrimaryTh << endl;
		nWinnersVec.clear();
		neAvalVec.clear();
	}
	/*
	while (hAmplification->GetMaximum() < 20) hAmplification->Rebin(2);
	hAmplification->Scale(1/hAmplification->GetMaximum());
	hAmplification->SetMaximum(1.2);
	hAmplification->SetMinimum(0);
	hAmplification->SetLineColor(kRed);
	 */
	TF1* fAmplification = GetFitCurve(hAmplification, false);
	fAmplification->SetLineColor(kRed);
	//hAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
	//hAmplification->GetXaxis()->SetRangeUser(0, 200000);
	
	hAmplification->Draw("hist same");
	fAmplification->Draw("same");
	hAmplification->SaveAs("rootFiles/Ar-iC4H10/model10/test.root");
	
	TLegend* lgd = new TLegend(0.3, 0.8, 0.9, 0.9);
	lgd->AddEntry(hFeAmplification, "Convolution with Fe spectrum", "l");
	lgd->AddEntry(hAmplification, "Fe source simulated", "l");
	lgd->Draw();
	
	cv->SaveAs(outputName);
	
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	
	return 0;
}


