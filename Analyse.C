#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>

#include "Include/Functions.C"
#include "Include/Transparency.C"


int Analyse(int modelNum, std::string gasName, std::vector<int> hvList) {
	
	bool extended = false;
	
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
	TString fSignalName = path + "signal";
	TString outputName = Form("Figures/model%d/analysis-%s", modelNum, gasName.c_str());
	TString outputName2 = Form("Figures/model%d/ionAnalysis-%s", modelNum, gasName.c_str());
	for (int k = 0; k< (int)hvList.size(); k++) {
		fileName += Form("-%d", hvList[k]);
		fSignalName += Form("-%d", hvList[k]);
		outputName += Form("-%d", hvList[k]);
		outputName2 += Form("-%d", hvList[k]);
	}
	fileName += ".root";
	fSignalName += ".root";
	outputName+=".pdf";
	outputName2+=".pdf";
	//outputName+=".root";
	
	std::cout << fSignalName << std::endl;
	
	std::map <std::string, int, NoSorting> electrode;
	LoadElectrodeMap(modelNum, electrode);
	
	int readoutElectrode = 0;
	int driftElectrode = 0;
	
	std::map<std::string, int>::iterator it = electrode.begin();
	for (it=electrode.begin(); it!=electrode.end(); ++it) {
		//std::cout << it->first << " => " << it->second << '\n';
		if (it->first == "pad") readoutElectrode = it->second;	// could be mesh, it depends on where you want to read
		else if (it->first == "drift") driftElectrode = it->second;
	}
	if (readoutElectrode == 0 || driftElectrode == 0) {std::cout << "Did not find drift or pad electrode" << std::endl; return 0;}
	
	/*
	 std::vector <Int_t> hv = {};
	 TObjArray* matches = TPRegexp( "-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?\\.root$" ).MatchS( fileName );
	 for (int k = 0; k<matches->GetLast(); k++) {
	 hv.push_back( (static_cast<TObjString*>(matches->At(k+1)))->GetString().Atoi() );
	 std::cout << hv[k] << std::endl;
	 }
	 */
	
	TCanvas* cv = new TCanvas("cv","cv", 1200, 1000);
	//cv->Divide(2);
	cv->Divide(3, 2);
	
	// Start with drawing the transparency
	cv->cd(1);
	DrawTransparency(modelNum, fSignalName);
	//return 0;
	
	// now draw the gain with the 2 different ways
	
	TFile* fSignal = TFile::Open(fSignalName, "READ");
	TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
	Int_t nAvalanche = tAvalanche->GetEntries();
	//Int_t nAmplification;
	vector <Int_t> *nWinnersVec = nullptr;
	vector <Int_t> gainVec = {};
	Int_t ni = 0, ionBackNum = 0;
	//tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
	tAvalanche->SetBranchAddress("amplificationElectrons", &nWinnersVec);
	tAvalanche->SetBranchAddress("ionNum", &ni);
	tAvalanche->SetBranchAddress("ionBackNum", &ionBackNum);
	
	// Create histogram of amplification electrons
	Int_t elMax = tAvalanche->GetMaximum("amplificationElectrons");
	TH1F* hAmplification = new TH1F("hAmplification", "Number of amplification electrons", elMax, 0, elMax);
	for (int k = 0; k< tAvalanche->GetEntries(); k++) {
		tAvalanche->GetEntry(k);
		gainVec = *nWinnersVec;
		//if (nAmplification>1) hAmplification->Fill(nAmplification);
		if (gainVec[0]>1) hAmplification->Fill(gainVec[0]);
	}
	while (hAmplification->GetMaximum() < 20) hAmplification->Rebin(2);
	hAmplification->Scale(1/hAmplification->GetMaximum());
	//hAmplification->SetMaximum(1.35);
	hAmplification->SetMaximum(1.15);
	hAmplification->SetLineColor(kBlue);
	TF1* fAmplification = GetFitCurve(hAmplification, false);
	fAmplification->SetLineColor(kBlue);
	
	Int_t iBinMax = hAmplification->GetMaximumBin();
	Double_t xMax = hAmplification->GetXaxis()->GetBinCenter( iBinMax );
	hAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hAmplification->GetRMS());
	
	/*
	 // Ignore convolution in the end
	 TFile* fConvoluted = TFile::Open(fileName, "READ");
	 TH1F* hFeAmplification = (TH1F*)fConvoluted->Get("hFeAmplification");
	 TH1F* hFeCharge = (TH1F*)fConvoluted->Get(Form("hFeCharge_%d", readoutElectrode));
	 while (hFeAmplification->GetMaximum() < 100) hFeAmplification->Rebin(2);
	 while (hFeCharge->GetMaximum() < 100) hFeCharge->Rebin(2);
	 //hFeAmplification->Rebin(8);
	 
	 hFeAmplification->Scale(1/hFeAmplification->GetMaximum());
	 hFeAmplification->SetMaximum(1.35);
	 //hFeAmplification->GetXaxis()->SetRangeUser(2, 10000);
	 hFeAmplification->SetLineColor(kBlue);
	 hFeAmplification->SetTitle("Gain with Fe source");
	 
	 TF1* f = GetFitCurve(hFeAmplification);
	 f->SetLineColor(kBlue);
	 
	 Int_t iBinMax = hFeAmplification->GetMaximumBin();
	 Double_t xMax = hFeAmplification->GetXaxis()->GetBinCenter( iBinMax );
	 hFeAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
	 
	 hFeCharge->Scale(1/hFeCharge->GetMaximum());
	 hFeCharge->SetLineColor(kRed);
	 TF1* f2 = GetFitCurve(hFeCharge);
	 f2->SetLineColor(kRed);
	 */
	
	
	// Now draw both spectra
	cv->cd(2);
	
	/*
	 // This is convolution
	 hFeAmplification->Draw("hist");
	 f->Draw("same");
	 // Add text to frame
	 TString txt = Form("Number of electrons --> Gain = %.0f #pm %.3f", f->GetParameter(0), f->GetParError(0));
	 
	 hFeCharge->Draw("hist same");
	 f2->Draw("same");
	 TString txt2 = Form("Induced charge --> Gain = %.0f #pm %.3f", f2->GetParameter(0), f2->GetParError(0));
	 
	 TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
	 legend->AddEntry(f,txt,"l");
	 legend->AddEntry(f2,txt2,"l");
	 legend->SetTextSize(0.04);
	 legend->Draw("same");
	 */
	
	hAmplification->Draw("hist");
	fAmplification->Draw("same");
	// Add text to frame
	TString txt = Form("Number of electrons --> Gain = %.0f #pm %.3f", fAmplification->GetParameter(0), fAmplification->GetParError(0));
	TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
	legend->AddEntry(fAmplification,txt,"l");
	legend->SetTextSize(0.04);
	//legend->Draw("same");
	
	// Write gain, sigma and res = sigma/gain on the plot
	double xPos = fAmplification->GetParameter(0) + 2.5*fAmplification->GetParameter(1);
	TLatex* txtGain = new TLatex(xPos, 0.9, Form("Gain = %.1f #pm %.1f ", fAmplification->GetParameter(0), fAmplification->GetParError(0)));
	TLatex* txtSigma = new TLatex(xPos, 0.8, Form("Sigma = %.1f #pm %.1f ", fAmplification->GetParameter(1), fAmplification->GetParError(1)));
	Double_t resolution = fAmplification->GetParameter(1)/fAmplification->GetParameter(0);
	Double_t resolutionError = resolution * TMath::Sqrt( Square(fAmplification->GetParError(0)/fAmplification->GetParameter(0)) + Square(fAmplification->GetParError(1)/fAmplification->GetParameter(1)) );
	std::string percent = "%";
	TLatex* txtRes = new TLatex(xPos, 0.7, Form("Sigma/Mean = %.1f #pm %.1f %s", resolution*100, resolutionError*100, percent.c_str()));
	
	txtGain->Draw("same"); txtSigma->Draw("same"); txtRes->Draw("same");
	
	
	std::cout << "\n\nStarting to draw the IBF now\n\n" << std::endl;
	
	// And finally draw the IBF
	
	// 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
	// 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
	// 3e étape : récupérer les ibf convolués
	// 4e étape : dessiner la courbe ibf = f(field ratio)
	
	// 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
	
	TH1F* hIbf = new TH1F("ibf", "ibf", 2000, 0, 100);
	for (int k = 0; k<nAvalanche; k++) {
		tAvalanche->GetEntry(k);
		tAvalanche->GetEntry(k);
		gainVec = *nWinnersVec;
		//if (nAmplification>1) {hIbf->Fill((double)ionBackNum/ni*100.);}
		if (gainVec[0]>1) {hIbf->Fill((double)ionBackNum/ni*100.);}
	}
	
	hIbf->SetXTitle("IBF (%)");
	hIbf->Scale(1/hIbf->GetMaximum());
	hIbf->SetMaximum(1.35);
	TF1* fIbf = GetFitCurve(hIbf);
	
	hIbf->GetXaxis()->SetRangeUser(0, 10.);
	if (fIbf->GetParameter(0) < 0.9) hIbf->GetXaxis()->SetRangeUser(0, 3.);
	if (fIbf->GetParameter(0) < 0.5) hIbf->GetXaxis()->SetRangeUser(0, 1.);
	bool gaussian = (fIbf->GetParameter(0)) > 0;
	if (!gaussian) fIbf = GetFitCurve(hIbf,false);
	
	// 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
	
	TTree* tChargeReadout = (TTree*)fSignal->Get(Form("tInducedCharge_%d", readoutElectrode));
	TTree* tChargeDrift = (TTree*)fSignal->Get(Form("tInducedCharge_%d", driftElectrode));
	Double_t chargeDrift, chargeReadout, ionChargeDrift, ionChargeReadout;
	tChargeReadout->SetBranchAddress("totalInducedCharge", &chargeReadout);
	tChargeDrift->SetBranchAddress("totalInducedCharge", &chargeDrift);
	tChargeReadout->SetBranchAddress("ionInducedCharge", &ionChargeReadout);
	tChargeDrift->SetBranchAddress("ionInducedCharge", &ionChargeDrift);
	
	int nChargeReadout = tChargeReadout->GetEntries();
	int nChargeDrift = tChargeDrift->GetEntries();
	if (nChargeReadout != nChargeDrift) {
		std::cout << "nChargeReadout != nChargeDrift" << std::endl;
		return 0;
	}
	int nCharge = nChargeDrift;
	
	int nBins = 10000;
	if (modelNum == 1 || modelNum == 16 || modelNum == 17 || modelNum == 18) nBins = 2000;
	TH1F* hIbfCharge = new TH1F("hIbfCharge", "hIbfCharge", nBins, 0, 100);
	TH1F* hIbfIonCharge = new TH1F("hIbfIonCharge", "hIbfIonCharge", nBins, 0, 100);
	for (int l = 0; l<nCharge; l++) {
		tChargeReadout->GetEntry(l);
		tChargeDrift->GetEntry(l);
		//tAvalanche->GetEntry(l);
		//if (nAmplification>1) {
		if (abs(chargeReadout)>10) {	// there is signal! it won't be 0 divided by 0
										// + if there's no signal, not relevant to compute an IBF
			hIbfCharge->Fill(abs(chargeDrift/chargeReadout*100.));
			hIbfIonCharge->Fill(abs(ionChargeDrift/ionChargeReadout*100.));
		}
	}
	hIbfCharge->Scale(1./hIbfCharge->GetMaximum());
	hIbfIonCharge->Scale(1./hIbfIonCharge->GetMaximum());
	TF1* fIbfCharge = GetFitCurve(hIbfCharge, gaussian);
	TF1* fIbfIonCharge = GetFitCurve(hIbfIonCharge, gaussian);
	
	/*
	 // 3e étape (nouvelle étape intermédiaire) : récupérer les histos IBF convolués, et les fitter
	 TH1F* hFeIbf = (TH1F*)fConvoluted->Get("hFeIbf");
	 TH1F* hFeIbfTotalCharge = (TH1F*)fConvoluted->Get("hFeIbfTotalCharge");
	 TH1F* hFeIbfIonCharge = (TH1F*)fConvoluted->Get("hFeIbfIonCharge");
	 hFeIbf->Scale(1/hFeIbf->GetMaximum());
	 hFeIbfTotalCharge->Scale(1/hFeIbfTotalCharge->GetMaximum());
	 hFeIbfIonCharge->Scale(1/hFeIbfIonCharge->GetMaximum());
	 
	 TF1* fFeIbf = GetFitCurve(hFeIbf);
	 TF1* fFeIbfTotalCharge = GetFitCurve(hFeIbfTotalCharge);
	 TF1* fFeIbfIonCharge = GetFitCurve(hFeIbfIonCharge);
	 hFeIbf->GetXaxis()->SetRangeUser(0, 10.);
	 */
	
	
	// 4e étape : dessiner les histos d'ibf
	cv->cd(3);
	hIbf->SetTitle("IBF");
	fIbf->SetLineColor(kBlue);
	hIbf->Draw("hist");
	fIbf->SetLineColor(kBlue);
	fIbf->Draw("same");
	
	/*
	 // Draw convoluted IBF
	 hFeIbf->SetLineColor(12);
	 fFeIbf->SetLineColor(12);
	 hFeIbf->Draw("hist same");
	 fFeIbf->Draw("same");
	 hFeIbfTotalCharge->SetLineColor(9);
	 fFeIbfTotalCharge->SetLineColor(9);
	 hFeIbfTotalCharge->Draw("hist same");
	 fFeIbfTotalCharge->Draw("same");
	 hFeIbfIonCharge->SetLineColor(8);
	 fFeIbfIonCharge->SetLineColor(8);
	 hFeIbfIonCharge->Draw("hist same");
	 fFeIbfIonCharge->Draw("same");
	 */
	
	// Draw not convoluted IBF histo with induced charge
	/*
	 hIbfCharge->SetLineColor(7);
	 fIbfCharge->SetLineColor(7);
	 hIbfCharge->Draw("hist same");
	 fIbfCharge->Draw("same");
	 */
	
	hIbfIonCharge->SetLineColor(6);
	fIbfIonCharge->SetLineColor(6);
	hIbfIonCharge->Draw("hist same");
	fIbfIonCharge->Draw("same");
	
	TString txtIbf1 = Form("IBF = %.3g #pm %.3f %s", fIbf->GetParameter(0), fIbf->GetParError(0), "%");
	TString txtIbf2 = Form("IBF = %.3g #pm %.3f %s", fIbfCharge->GetParameter(0), fIbfCharge->GetParError(0), "%");
	TString txtIbf3 = Form("IBF = %.3g #pm %.3f %s", fIbfIonCharge->GetParameter(0), fIbfIonCharge->GetParError(0), "%");
	TLegend* legend2 = new TLegend(0.1,0.75,0.9,0.9);
	legend2->AddEntry(hIbf,"IBF ratio --> " + txtIbf1,"l");
	//legend2->AddEntry(hIbfCharge,"All induced charges --> " + txtIbf2,"l");
	legend2->AddEntry(hIbfIonCharge,"Induced ion charges --> " + txtIbf3,"l");
	legend2->SetTextSize(0.04);
	legend2->Draw("same");
	
	
	/*
	 cv->cd(4);
	 hIbfCharge->SetTitle("IBF (zoom)");
	 hIbfCharge->SetXTitle("IBF (%)");
	 hIbfCharge->GetXaxis()->SetRangeUser(0,0.4);
	 hIbfCharge->SetMaximum(1.3);
	 hIbfCharge->Draw("hist");
	 fIbfCharge->Draw("same");
	 // Draw convoluted IBF
	 hFeIbf->Draw("hist same");
	 hFeIbfTotalCharge->Draw("hist same");
	 hFeIbfIonCharge->Draw("hist same");
	 
	 hIbfIonCharge->Draw("hist same");
	 fIbfIonCharge->Draw("same");
	 
	 TLegend* legend3 = (TLegend*)legend2->Clone();
	 legend3->DeleteEntry();
	 legend3->Draw();
	 */
	
	//cv->cd(5);
	cv->cd(4);
	DrawDetector(modelNum, hvList);
	
	TText* txtGas = new TText(.4,.95,Form("Gas: %s", gasName.c_str()));
	txtGas->Draw("same");
	
	// Draw distribution of where electrons are created
	cv->cd(5);
	std::vector<float> *electronStartPointsInput = 0;
	tAvalanche->SetBranchAddress("electronStartPoints", &electronStartPointsInput);
	
	std::vector<float> electronStartPoints = {};
	TH1F* zElDistribution = new TH1F("hZelectrons", "Start z of electrons", 1000, 0, damp*1.2);
	for (int k = 0; k < nAvalanche; k++) {
		tAvalanche->GetEntry(k);
		electronStartPoints = *electronStartPointsInput;
		for (int j = 0; j< (int)electronStartPoints.size(); j++) {
			zElDistribution->Fill(electronStartPoints[j]);
		}
		electronStartPoints.clear();
	}
	zElDistribution->Scale(1/zElDistribution->GetMaximum());
	if (!(modelNum==1 || (modelNum >15 && modelNum < 19)) ) zElDistribution->SetMaximum(0.05);
	zElDistribution->GetXaxis()->SetTitle("z (cm)");
	zElDistribution->Draw("hist");
	
	cv->cd(6);
	DrawDyingIons(modelNum, fSignalName);
	
	cv->SaveAs(outputName);
	
	if (!extended) {
		time_t t1 = time(NULL);
		PrintTime(t0, t1);
		
		return 0;
	}
	
	// Draw histograms of IBF (see in Transparency.C)
	
	TCanvas* cv2 = new TCanvas("cv2","cv2", 1200, 1000);
	//cv->Divide(2);
	cv2->Divide(3, 2);
	
	cv2->cd(1);
	DrawNionsInDriftRegion(fSignalName);
	
	cv2->cd(2);
	DrawDyingIons(modelNum, fSignalName);
	cv2->cd(3);
	DrawDyingIons2d(modelNum, fSignalName);
	cv2->cd(4);
	DrawDyingIons3d(modelNum, fSignalName);
	
	cv2->SaveAs(outputName2);
	
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	
	return 0;
}


