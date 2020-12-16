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

#include "Functions.C"
#include "Utils.C"
#include "Geometry.C"
#include "Transparency.C"

using namespace std;

void DrawAmplificationElectrons(string gasName = "Ar-iC4H10", TString fSignalName="", bool useFeSource = false) {
	/* Draw the gain */
	
	TFile* fSignal = TFile::Open(fSignalName, "READ");
	TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
	Int_t nAvalanche = tAvalanche->GetEntries();
	
	if (useFeSource) {
		cout << "There are " << (int)tAvalanche->GetEntries() << " entries in this file" << endl;
	}
	
	Int_t nPrimaryTh = GetPrimary(gasName);
	
	// Create histogram of amplification electrons
	Int_t elMax = GetMaxAmp(*tAvalanche);	// must be called before setting the addresses to the tree
	if (useFeSource) elMax = int((double)elMax/nPrimaryTh);
	
	/*
	 Int_t nAmplification;
	 tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
	 */
	vector<Int_t> *neAvalVecIn = nullptr, *nWinnersVecIn = nullptr;
	tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVecIn);
	tAvalanche->SetBranchAddress("avalancheSize", &neAvalVecIn);
	
	vector<Int_t> nWinnersVec = {}, neAvalVec = {};
	TH1F* hAmplification = new TH1F("hAmplification", "Number of amplification electrons", elMax, 0, 2*elMax);
	const int nEntries = (int)tAvalanche->GetEntries();
	for (int k = 0; k< (int)tAvalanche->GetEntries(); k++) {
		tAvalanche->GetEntry(k);
		nWinnersVec = *nWinnersVecIn;
		neAvalVec = *neAvalVecIn;
		int nAmplification = 0;
		//cout << nWinnersVec.size() << endl;
		for (int l = 0; l< (int)nWinnersVec.size(); l++) {nAmplification += nWinnersVec[l];}
		//cout << nAmplification << endl;
		if (useFeSource) nAmplification /= (double)nPrimaryTh;
		if (nAmplification > 5) hAmplification->Fill(nAmplification);
		//cout << (double)nAmplification/nPrimaryTh << endl;
		nWinnersVec.clear();
		neAvalVec.clear();
	}
	while (hAmplification->GetMaximum() < 30) hAmplification->Rebin(2);
	if ( nEntries > 1000) {while (hAmplification->GetMaximum() < 70) hAmplification->Rebin(2);}
	if ( nEntries > 5000) {while (hAmplification->GetMaximum() < 120) hAmplification->Rebin(2);}
	if ( nEntries > 10000) {while (hAmplification->GetMaximum() < 200) hAmplification->Rebin(2);}
	hAmplification->Scale(1/hAmplification->GetMaximum());
	hAmplification->SetMaximum(1.3);
	hAmplification->SetMinimum(0);
	hAmplification->SetLineColor(kBlue);
	
	TF1* fAmplification = GetFitCurve(hAmplification, useFeSource);	// if Fe source, fit = gauss; else fit = landau
	fAmplification->SetLineColor(kBlue);
	
	Int_t iBinMax = hAmplification->GetMaximumBin();
	Double_t xMax = hAmplification->GetXaxis()->GetBinCenter( iBinMax );
	hAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hAmplification->GetRMS());
	
	if (xMax > 8000) hAmplification->GetXaxis()->SetMaxDigits(3);
	
	hAmplification->Draw("hist");
	fAmplification->Draw("same");
	
	// Add text to frame
	TString txt = Form("Number of electrons --> Gain = %.0f #pm %.3f", fAmplification->GetParameter(0), fAmplification->GetParError(0));
	TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
	legend->AddEntry(fAmplification,txt,"l");
	legend->SetTextSize(0.04);
	//legend->Draw("same");
	
	// Write gain, sigma and res = sigma/gain on the plot
	double xPos = fAmplification->GetParameter(0)*0.2;
	TLatex* txtGain = new TLatex(xPos, 1.1, Form("Gain = %.1f #pm %.1f ", fAmplification->GetParameter(0), fAmplification->GetParError(0)));
	TLatex* txtSigma = new TLatex(xPos, 1.0, Form("Sigma = %.1f #pm %.1f ", fAmplification->GetParameter(1), fAmplification->GetParError(1)));
	Double_t resolution = fAmplification->GetParameter(1)/fAmplification->GetParameter(0);
	Double_t resolutionError = resolution * TMath::Sqrt( Square(fAmplification->GetParError(0)/fAmplification->GetParameter(0)) + Square(fAmplification->GetParError(1)/fAmplification->GetParameter(1)) );
	std::string percent = "%";
	TLatex* txtRes = new TLatex(xPos, 0.9, Form("Sigma/Mean = %.1f #pm %.1f %s", resolution*100, resolutionError*100, percent.c_str()));
	
	txtGain->Draw("same"); txtSigma->Draw("same"); txtRes->Draw("same");
}


void DrawFeConvolution(TString fConvolutedName="") {
	// First draw convoluted spectrum
	TFile* fConvoluted = TFile::Open(fConvolutedName, "READ");
	
	TH1F* hFeAmplification = (TH1F*)fConvoluted->Get("hFeAmplification");
	while (hFeAmplification->GetMaximum() < 100) hFeAmplification->Rebin(2);
	
	hFeAmplification->Scale(1/hFeAmplification->GetMaximum());
	hFeAmplification->SetMaximum(1.3);
	//hFeAmplification->GetXaxis()->SetRangeUser(2, 10000);
	hFeAmplification->SetLineColor(kBlue);
	hFeAmplification->SetTitle("Gain convoluted with Fe source");
	
	TF1* f = GetFitCurve(hFeAmplification);
	f->SetLineColor(kBlue);
	
	Int_t iBinMax = hFeAmplification->GetMaximumBin();
	Double_t xMax = hFeAmplification->GetXaxis()->GetBinCenter( iBinMax );
	hFeAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
	
	if (xMax > 8000) hFeAmplification->GetXaxis()->SetMaxDigits(3);
	
	hFeAmplification->Draw("hist");
	f->Draw("same");
	
	// Write gain, sigma and res = sigma/gain on the plot
	double xPos = f->GetParameter(0)*0.2;
	TLatex* txtGain = new TLatex(xPos, 1.1, Form("Gain = %.1f #pm %.1f ", f->GetParameter(0), f->GetParError(0)));
	TLatex* txtSigma = new TLatex(xPos, 1.0, Form("Sigma = %.1f #pm %.1f ", f->GetParameter(1), f->GetParError(1)));
	Double_t resolution = f->GetParameter(1)/f->GetParameter(0);
	Double_t resolutionError = resolution * TMath::Sqrt( Square(f->GetParError(0)/f->GetParameter(0)) + Square(f->GetParError(1)/f->GetParameter(1)) );
	std::string percent = "%";
	TLatex* txtRes = new TLatex(xPos, 0.9, Form("Sigma/Mean = %.1f #pm %.1f %s", resolution*100, resolutionError*100, percent.c_str()));
	
	txtGain->Draw("same"); txtSigma->Draw("same"); txtRes->Draw("same");
	
	/*
	 TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
	 legend->AddEntry(f,txt,"l");
	 legend->AddEntry(f2,txt2,"l");
	 legend->SetTextSize(0.04);
	 legend->Draw("same");
	 */
}


void DrawFeChargeConvolution(TString fConvolutedName="", int readoutElectrode = 0) {
	
	TFile* fConvoluted = TFile::Open(fConvolutedName, "READ");
	TH1F* hFeCharge = (TH1F*)fConvoluted->Get(Form("hFeCharge_%d", readoutElectrode));
	while (hFeCharge->GetMaximum() < 100) hFeCharge->Rebin(2);
	//hFeAmplification->Rebin(8);
	
	hFeCharge->Scale(1/hFeCharge->GetMaximum());
	hFeCharge->SetLineColor(kRed);
	TF1* f2 = GetFitCurve(hFeCharge);
	f2->SetLineColor(kRed);
	
	hFeCharge->Draw("hist same");
	f2->Draw("same");
	TString txt2 = Form("Induced charge --> Gain = %.0f #pm %.3f", f2->GetParameter(0), f2->GetParError(0));
	/*
	 TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
	 legend->AddEntry(f,txt,"l");
	 legend->AddEntry(f2,txt2,"l");
	 legend->SetTextSize(0.04);
	 legend->Draw("same");
	 */
}

int Spectrum(int modelNum, std::string gasName, std::vector<int> hvList, bool useFeSource) {
	
	time_t t0 = time(NULL);
	gStyle->SetTitleFontSize(.06);
	gStyle->SetTitleSize(.06);
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetLabelSize(.04, "XY");
	gStyle->SetMarkerSize(0.3);
	gStyle->SetTextSize(0.05);
	
	int electrodeNum = GetElectrodeNum(modelNum);
	if ((int)hvList.size() != electrodeNum-1) {cout << "Wrong hv input" << endl; return 0;}
	
	//Load geometry parameters
	double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
	
	// input and output files
	const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	TString fConvolutedName = path + "fe-spectrum-convoluted";
	TString fSignalName = path;
	if (useFeSource) fSignalName += "fe-signal-noibf";
	else fSignalName += "signal-noibf";
	//TString fSignalName = path + "signal";
	TString outputName = Form("Figures/model%d/spectra-%s", modelNum, gasName.c_str());
	for (int k = 0; k< (int)hvList.size(); k++) {
		fConvolutedName += Form("-%d", hvList[k]);
		fSignalName += Form("-%d", hvList[k]);
		outputName += Form("-%d", hvList[k]);
	}
	fConvolutedName += ".root";
	fSignalName += ".root";
	outputName+=".pdf";
	
	cout << fSignalName << endl;
	if (useFeSource) {cout << endl << "USING FE SOURCE" << endl << endl;}
	
	
	TCanvas* cv = new TCanvas("cv","cv", 1200, 1000);
	cv->Divide(2);
	
	cv->cd(1);
	DrawAmplificationElectrons(gasName, fSignalName, useFeSource);
	
	cv->cd(2);
	DrawFeConvolution(fConvolutedName);
	
	/*
	 TLegend* lgd = new TLegend(0.3, 0.8, 0.9, 0.9);
	 lgd->AddEntry(hFeAmplification, "Convolution with Fe spectrum", "l");
	 //lgd->AddEntry(hAmplification, "Fe source simulated", "l");
	 lgd->Draw();
	 */
	
	cv->SaveAs(outputName);
	
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	
	return 0;
}


