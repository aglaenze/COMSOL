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


Double_t FitFunctionCrystalBall( Double_t* x, Double_t* par ) { //(x, alpha, n sigma, mu)
	return par[0]*ROOT::Math::crystalball_pdf( x[0], par[1], par[2], par[3], par[4] ); }

Double_t FitGauss( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*TMath::Gaus( x[0], par[0], par[1]); }

//____________________________________________
Double_t FitFunctionExp( Double_t* x, Double_t* par ) {
	return TMath::Exp( par[0] + par[1]*x[0] ); }

Double_t FitLandau( Double_t* x, Double_t* par ) { //Double_t TMath::Landau	(Double_t x, Double_t mu = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*TMath::Landau( x[0], par[0], par[1]); }

TF1* GetFitGain(TH1F* h) {
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	std::cout << "xMax = " << xMax << std::endl;
	std::cout << "maximum = " << h->GetMaximum() << std::endl;
	
	Int_t fitRangeMin = xMax - 1.1 * h->GetRMS();
	Int_t fitRangeMax = xMax + 1.1 * h->GetRMS();
	
	TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("Mean", "Sigma", "Amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}

TF1* GetFitIbf(TH1F* h, bool gauss = true) {
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	std::cout << "xMax = " << xMax << std::endl;
	std::cout << "maximum = " << h->GetMaximum() << std::endl;
	
	Int_t fitRangeMin = 0;
	Int_t fitRangeMax = xMax + h->GetRMS();
	if (gauss) {
		fitRangeMin = xMax - h->GetRMS();
		fitRangeMax = xMax + 4*h->GetRMS();
	}
	TF1* f;
	if (gauss) f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
	else f = new TF1( "FitFunction", FitLandau, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("Mean", "Sigma", "Amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	std::cout << "\n\nh->GetRMS() = " << h->GetRMS() << std::endl;
	std::cout << "\n\nh->GetMaximum() = " << h->GetMaximum() << std::endl;
	if (!gauss) {
		f->SetParLimits(0, 0.8*xMax, 1.3*xMax);
		f->SetParLimits(1,0, 0.001*h->GetRMS());
		f->FixParameter(2, 1);
		//f->SetParLimits(2, 0.5*h->GetMaximum(), 5*h->GetMaximum());
		//f->SetParameter(2,4*h->GetMaximum());
	}
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}


int Analyse(int modelNum, std::string gasName, std::vector<int> hvList) {
	

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
	if (hvList.size() != electrodeNum-1) {std::cout << "Wrong hv input" << std::endl; return 0;}
	
	// input and output files
	const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	TString fileName = path + "fe-spectrum-convoluted";
	TString fSignalName = path + "signal";
	TString outputName = Form("Figures/model%d/analysis-%s", modelNum, gasName.c_str());
	for (int k = 0; k< hvList.size(); k++) {
		fileName += Form("-%d", hvList[k]);
		fSignalName += Form("-%d", hvList[k]);
		outputName += Form("-%d", hvList[k]);
	}
	fileName += ".root";
	fSignalName += ".root";
	outputName+=".pdf";
	
	std::cout << fSignalName << std::endl;
	
	std::map <std::string, int> electrode;
	LoadElectrodeMap(modelNum, electrode);
	
	int readoutElectrode = electrode["pad"]; // could be mesh, it depends on where you want to read
	int driftElectrode = electrode["drift"];
	
	/*
	std::vector <Int_t> hv = {};
	TObjArray* matches = TPRegexp( "-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?\\.root$" ).MatchS( fileName );
	for (int k = 0; k<matches->GetLast(); k++) {
		hv.push_back( (static_cast<TObjString*>(matches->At(k+1)))->GetString().Atoi() );
		std::cout << hv[k] << std::endl;
	}
	 */
	
	TCanvas* cv = new TCanvas("cv","cv", 800, 1000);
	//cv->Divide(2);
	cv->Divide(2, 2);
	
	// Start with drawing the transparency
	
	TFile* fSignal = TFile::Open(fSignalName, "READ");
	TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
	Int_t nAvalanche = tAvalanche->GetEntries();
	
	Int_t nAmplification;
	Int_t ni = 0, ionBackNum = 0;
	tAvalanche->SetBranchAddress("amplificationElectrons", &nAmplification);
	tAvalanche->SetBranchAddress("ionNum", &ni);
	tAvalanche->SetBranchAddress("ionBackNum", &ionBackNum);
	
	Int_t winners = 0;
	TH1F* hTransparency = new TH1F("hTransparency", "hTransparency", 4, -1, 3);
	
	// On récupère en même temps les fichiers d'IBF et de transparence
	// 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
	
	TH1F* hIbf = new TH1F("ibf", "ibf", 2000, 0, 100);
	for (int k = 0; k<nAvalanche; k++) {
		tAvalanche->GetEntry(k);
		if (nAmplification>1) {
			hIbf->Fill((double)ionBackNum/ni*100.);
			hTransparency->Fill(1);
			winners++;
		}
		else {hTransparency->Fill(0);}
		//if (ibfRatio*100. > 100.) std::cout << "problem !" << std::endl;
	}
	
	cv->cd(1);
	hTransparency->SetTitle("Transparency");
	hTransparency->Scale(1./nAvalanche);
	hTransparency->SetMaximum(1.2);
	hTransparency->Draw("hist");
	// Write text with the value of transparency
	Double_t transp = (double)winners/nAvalanche;
	
	TText* txttr = new TText(0.2,0.9*hTransparency->GetMaximum(),Form("Transparency = %.1f %s", transp*100, "%"));
	txttr->Draw();
	
	// now draw the gain with the 2 different ways
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
	
	TF1* f = GetFitGain(hFeAmplification);
	f->SetLineColor(kBlue);
	
	Int_t iBinMax = hFeAmplification->GetMaximumBin();
	Double_t xMax = hFeAmplification->GetXaxis()->GetBinCenter( iBinMax );
	hFeAmplification->GetXaxis()->SetRangeUser(0, xMax + 3*hFeAmplification->GetRMS());
	
	hFeCharge->Scale(1/hFeCharge->GetMaximum());
	hFeCharge->SetLineColor(kRed);
	TF1* f2 = GetFitGain(hFeCharge);
	f2->SetLineColor(kRed);
	
	
	// Now draw both spectra
	cv->cd(2);
	
	hFeAmplification->Draw("hist");
	f->Draw("same");
	// Add text to frame
	TString txt = Form("Number of electrons --> Gain = %.0f #pm %.3f", f->GetParameter(0), f->GetParError(0));
	//Latex* txt = new TLatex(0.5*xMax,1,Form("Number of electrons --> Gain = %.3g #pm %.3f", f->GetParameter(0), f->GetParError(0)));
	
	hFeCharge->Draw("hist same");
	f2->Draw("same");
	TString txt2 = Form("Induced charge --> Gain = %.0f #pm %.3f", f2->GetParameter(0), f2->GetParError(0));
	
	TLegend* legend = new TLegend(0.1,0.75,0.9,0.9);
	legend->AddEntry(f,txt,"l");
	legend->AddEntry(f2,txt2,"l");
	legend->SetTextSize(0.04);
	legend->Draw("same");
	
	std::cout << "\n\nStarting to draw the IBF now\n\n" << std::endl;
	
	// And finally draw the IBF
	
	// 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
	// 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
	// 3e étape : récupérer les ibf convolués
	// 4e étape : dessiner la courbe ibf = f(field ratio)
	
	hIbf->SetXTitle("IBF (%)");
	hIbf->Scale(1/hIbf->GetMaximum());
	hIbf->SetMaximum(1.35);
	TF1* fIbf = GetFitIbf(hIbf);
	
	Int_t iBinMax2 = hIbf->GetMaximumBin();
	Double_t xMax2 = hIbf->GetXaxis()->GetBinCenter( iBinMax2 );
	//hIbf->GetXaxis()->SetRangeUser(0, xMax2 + 3*hIbf->GetRMS());
	hIbf->GetXaxis()->SetRangeUser(0, 10.);

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
	hIbfCharge->Scale(1/hIbfCharge->GetMaximum());
	hIbfIonCharge->Scale(1/hIbfIonCharge->GetMaximum());
	TF1* fIbfCharge = GetFitIbf(hIbfCharge, false);
	TF1* fIbfIonCharge = GetFitIbf(hIbfIonCharge, false);
	
	// 3e étape (nouvelle étape intermédiaire) : récupérer les histos IBF convolués, et les fitter
	TH1F* hFeIbf = (TH1F*)fConvoluted->Get("hFeIbf");
	TH1F* hFeIbfTotalCharge = (TH1F*)fConvoluted->Get("hFeIbfTotalCharge");
	TH1F* hFeIbfIonCharge = (TH1F*)fConvoluted->Get("hFeIbfIonCharge");
	hFeIbf->Scale(1/hFeIbf->GetMaximum());
	hFeIbfTotalCharge->Scale(1/hFeIbfTotalCharge->GetMaximum());
	hFeIbfIonCharge->Scale(1/hFeIbfIonCharge->GetMaximum());
	
	TF1* fFeIbf = GetFitIbf(hFeIbf);
	TF1* fFeIbfTotalCharge = GetFitIbf(hFeIbfTotalCharge);
	TF1* fFeIbfIonCharge = GetFitIbf(hFeIbfIonCharge);
	hFeIbf->GetXaxis()->SetRangeUser(0, 10.);
	
	
	// 4e étape : dessiner les histos d'ibf
	cv->cd(3);
	hIbf->SetTitle("IBF");
	fIbf->SetLineColor(kBlue);
	hIbf->Draw("hist");
	fIbf->SetLineColor(kBlue);
	fIbf->Draw("same");

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
	
	// Draw not convoluted IBF histo with induced charge
	hIbfCharge->SetLineColor(7);
	fIbfCharge->SetLineColor(7);
	hIbfCharge->Draw("hist same");
	fIbfCharge->Draw("same");
	
	hIbfIonCharge->SetLineColor(6);
	fIbfIonCharge->SetLineColor(6);
	hIbfIonCharge->Draw("hist same");
	fIbfIonCharge->Draw("same");
	
	TString txtIbf1 = Form("IBF = %.3g #pm %.3f %s", fIbf->GetParameter(0), fIbf->GetParError(0), "%");
	TString txtIbf2 = Form("IBF = %.3g #pm %.3f %s", fIbfCharge->GetParameter(0), fIbfCharge->GetParError(0), "%");
	TString txtIbf3 = Form("IBF = %.3g #pm %.3f %s", fIbfIonCharge->GetParameter(0), fIbfIonCharge->GetParError(0), "%");
	TLegend* legend2 = new TLegend(0.1,0.7,0.9,0.9);
	legend2->AddEntry(hIbf,"IBF ratio --> " + txtIbf1,"l");
	legend2->AddEntry(hIbfCharge,"All induced charges --> " + txtIbf2,"l");
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
	
	TText* txtGas = new TText(.35,.7,Form("Gas: %s", gasName.c_str()));
	//txtGas->Draw();
	txtGas->Draw("same");
	

	cv->SaveAs(outputName);
	time_t t1 = time(NULL);
	PrintTime(t0, t1);
	
	return 0;
}


