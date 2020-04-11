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


Double_t FitGauss( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*TMath::Gaus( x[0], par[0], par[1]); }

Double_t FitLandau( Double_t* x, Double_t* par ) { //Double_t TMath::Landau	(Double_t x, Double_t mu = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*TMath::Landau( x[0], par[0], par[1]); }

TF1* GetFitIbf(TH1F* h, bool gauss = true) {
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	std::cout << "xMax = " << xMax << std::endl;
	std::cout << "maximum = " << h->GetMaximum() << std::endl;
	
	Int_t fitRangeMin = xMax - h->GetRMS();
	Int_t fitRangeMax = xMax + h->GetRMS();
	
	TF1* f;
	if (gauss) f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
	else f = new TF1( "FitFunction", FitLandau, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("Mean", "Sigma", "Amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	f->SetParLimits(0,0,10);
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}


int GetIbf() {
	
	//______________________
	// variables
	std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
									   //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
	const int modelNum = 1;
	//____________________
	
	time_t t0 = time(NULL);
	
	gStyle->SetTitleFontSize(.06);
	gStyle->SetTitleSize(.06);
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetLabelSize(.05, "XY");
	
	const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
	
	// Get number of files to look at
	Int_t num = 6;
	//Int_t num = GetNumberOfFiles(path, "signal");
	Int_t num2 = int(num/2.);
	if (num/2.> num2) num2+=1;
	TCanvas* c2 = new TCanvas("c2", "c2", 400, 300*num2);
	c2->Divide(2, num2);
	
	Double_t ibfList[num], ibfErrorList[num], ibfList2[num], ibfErrorList2[num], ibfList3[num];
	Double_t hvMeshList[num], hvDriftList[num], hvRatioList[num];
	double damp = 0.0128, ddrift = 0.5;
	
	const Int_t dataNum = dataQuantity(gasName);
	
	std::map <std::string, int> electrode;
	LoadElectrodeMap(modelNum, electrode);
	
	int readoutElectrode = electrode["pad"];    // could be mesh, it depends on where you want to read
	int driftElectrode = electrode["drift"];
	
	for (int k = 0; k<num; k++) {
		int hvMesh = 340+20*k;
		int hvDrift = 200+hvMesh;
		hvRatioList[k] = (double)hvMesh/(hvDrift-hvMesh)*(ddrift-damp)/damp;
		
		TFile* fSignal = TFile::Open(Form("rootFiles/Ar-iC4H10/model1/signal-%d-%d.root", hvMesh, hvDrift), "READ");
		
		// 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
		// 2e étape : récupérer les fichiers de charge induite, dessiner le ratio dans un TH1, et les fitter
		// 3e étape: faire une ligne verticale de l'estimation du ratio des intégrales
		// 4e étape : dessiner la courbe ibf = f(field ratio)
		
		// 1ere étape: récupérer les fichiers de ibf, les dessiner dans un TH1, et les fitter
		TTree* tAvalanche = (TTree*)fSignal->Get("tAvalanche");
		Double_t ibfRatio;
		tAvalanche->SetBranchAddress("ibfRatio", &ibfRatio);
		int nAvalanche = tAvalanche->GetEntries();
		
		TH1F* hIbf = new TH1F("ibf", "ibf", 1000, 0, 100);
		hIbf->SetXTitle("IBF (%)");
		for (int l = 0; l<nAvalanche; l++) {
			tAvalanche->GetEntry(l);
			hIbf->Fill(ibfRatio*100.);
		}
		
		hIbf->Scale(1/hIbf->GetMaximum());
		hIbf->SetMaximum(1.2);
		TF1* f = GetFitIbf(hIbf);
		
		
		Int_t iBinMax = hIbf->GetMaximumBin();
		Double_t xMax = hIbf->GetXaxis()->GetBinCenter( iBinMax );
		hIbf->GetXaxis()->SetRangeUser(0, xMax + 3*hIbf->GetRMS());
		
		hIbf->GetXaxis()->SetRangeUser(0, 10.);
		hvMeshList[k] = hvMesh;
		ibfList[k] = f->GetParameter(0);
		ibfErrorList[k] = f->GetParError(0);
		
		
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
		
		TH1F* hIbfCharge = new TH1F("hIbfCharge", "hIbfCharge", 1000, 0, 100);
		TH1F* hIbfIonCharge = new TH1F("hIbfIonCharge", "hIbfIonCharge", 1000, 0, 100);
		for (int l = 0; l<nCharge; l++) {
			tChargeReadout->GetEntry(l);
			tChargeDrift->GetEntry(l);
			hIbfCharge->Fill(abs(chargeDrift/chargeReadout*100.));
			hIbfIonCharge->Fill(abs(ionChargeDrift/ionChargeReadout*100.));
		}
		hIbfCharge->Scale(1/hIbfCharge->GetMaximum());
		hIbfCharge->SetMaximum(1.2);
		hIbfIonCharge->Scale(1/hIbfIonCharge->GetMaximum());
		TF1* f2 = GetFitIbf(hIbfCharge, false);
		TF1* f3 = GetFitIbf(hIbfIonCharge, false);
		
		
		/*
		 Int_t iBinMax2 = hIbfCharge->GetMaximumBin();
		 Double_t xMax2 = hIbfCharge->GetXaxis()->GetBinCenter( iBinMax );
		 hIbfCharge->GetXaxis()->SetRangeUser(0, xMax + 3*hIbf->GetRMS());
		 */
		
		hIbfCharge->GetXaxis()->SetRangeUser(0, 10.);
		ibfList2[k] = f2->GetParameter(0);
		ibfErrorList2[k] = f2->GetParError(0);
		
		// 3e étape: dessiner une ligne verticale de l'estimation du ratio des intégrales
		
		TTree* tSignalReadout = (TTree*)fSignal->Get(Form("tSignal_%d", readoutElectrode));
		TTree* tSignalDrift = (TTree*)fSignal->Get(Form("tSignal_%d", driftElectrode));
		Double_t currentDriftIntMax = tSignalDrift->GetMaximum("totalCurrentInt");
		Double_t currentReadoutIntMax = tSignalReadout->GetMaximum("totalCurrentInt");
		Double_t currentReadoutIntMin = tSignalReadout->GetMinimum("totalCurrentInt");
		if (currentReadoutIntMin<-100) currentReadoutIntMax = -tSignalReadout->GetMinimum("totalCurrentInt");
		
		Double_t ibfEstimation = 100*currentDriftIntMax/currentReadoutIntMax;
		ibfList3[k] = ibfEstimation;
		Double_t yMax = 10;
		TLine* ibfLine = new TLine(ibfEstimation, 0, ibfEstimation, yMax);
		ibfLine->SetLineColor(2);
		
		
		// 4e étape : dessiner les courbes ibf = f(field ratio)
		c2->cd(k+1);
		// Upper plot will be in pad1
		TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 0.95);
		pad1->SetTopMargin(0.2);
		pad1->SetBottomMargin(0); // Upper and lower plot are joined
		pad1->SetLeftMargin(0.15);
		//pad1->SetGridx();         // Vertical grid
		pad1->Draw();             // Draw the upper pad: pad1
		pad1->cd();               // pad1 becomes the current pad
		hIbf->Draw("hist");
		f->Draw("same");
		
		
		// Add text to frame
		TLatex* txt = new TLatex(1.2*xMax,1,Form("IBF = %.3g #pm %.3f %s", f->GetParameter(0), f->GetParError(0), "%"));
		TLatex* txt2 = new TLatex(1.2*xMax,0.9,Form("Field ratio = %.2f", hvRatioList[k]));
		txt->SetTextSize(0.05) ;
		txt2->SetTextSize(0.05) ;
		txt->Draw();
		txt2->Draw();
		
		// lower plot will be in pad
		c2->cd(k+1);          // Go back to the main canvas before defining pad2
		TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.5);
		pad2->SetTopMargin(0);
		pad2->SetLeftMargin(0.15);
		pad2->SetBottomMargin(0.2);
		//pad2->SetGridx(); // vertical grid
		pad2->Draw();
		pad2->cd();       // pad2 becomes the current pad
		hIbfCharge->Draw("hist");
		hIbfIonCharge->SetLineColor(6);
		hIbfIonCharge->Draw("hist same");
		f2->Draw("same");
		
		TLegend* legend = new TLegend(0.6,0.7,0.9,0.9);
		legend->AddEntry(hIbfCharge,"All induced charges","l");
		legend->AddEntry(hIbfIonCharge,"Induced ion charges","l");
		legend->Draw("same");
		
		TLatex* txt3 = new TLatex(1.2*xMax,1,Form("IBF = %.3g #pm %.3f %s", f2->GetParameter(0), f2->GetParError(0), "%"));
		TLatex* txt4 = new TLatex(1.2*xMax,0.9,Form("Field ratio = %.2f", hvRatioList[k]));
		TLatex* txtEstim = new TLatex(1.2*xMax,0.5,Form("IBF estimation = %.2f %s (current integration)", ibfEstimation, "%"));
		txt3->SetTextSize(0.05) ;
		txt4->SetTextSize(0.05) ;
		txt3->Draw();
		txt4->Draw();
		txtEstim->Draw();
		
		
		c2->cd(k+1);
		//ibfLine->Draw("same");
	}
	c2->SaveAs(Form("Figures/ibf_%s_model%d.pdf", gasName.c_str(), modelNum));
	
	
	TCanvas* c5 = new TCanvas("c5", "c5", 200, 300);
	c5->SetGrid();
	c5->Divide(1,2);
	
	c5->cd(1);
	
	TGraphErrors* grSim1 = new TGraphErrors(num, hvMeshList, ibfList, 0, ibfErrorList);
	grSim1->SetTitle("IBF curve in the Micromegas");
	grSim1->GetXaxis()->SetTitle( "V_{mesh}" );
	grSim1->GetYaxis()->SetTitle( "IBF (%)" );
	//grSim1->GetXaxis()->SetLimits(0, 4);  // along X axis
	grSim1->GetHistogram()->SetMinimum(0.);   // along Y axis
	grSim1->GetHistogram()->SetMaximum(4.);   // along Y axis
	grSim1->SetMarkerStyle(20);
	grSim1->SetMarkerSize(0.3);
	grSim1->SetMarkerColor(1);
	//grSim1->GetXaxis()->SetLimits(hvMeshList[0]-5, hvMeshList[num]+5);
	grSim1->Draw("ALP");
	
	TGraphErrors* grSim2 = new TGraphErrors(num, hvMeshList, ibfList2, 0, ibfErrorList2);
	grSim2->SetMarkerStyle(20);
	grSim2->SetMarkerSize(0.3);
	grSim2->SetMarkerColor(2);
	grSim2->Draw("LP same");
	
	TGraph* grSim3 = new TGraph(num, hvMeshList, ibfList3);
	grSim3->SetMarkerStyle(20);
	grSim3->SetMarkerSize(0.3);
	grSim3->SetMarkerColor(3);
	grSim3->Draw("LP same");
	
	// Same with data
	//const Int_t dataNum = dataQuantity(gasName);
	Double_t hvMeshListIBF1[dataNum], hvDriftListIBF1[dataNum], ionBackFlowCorrectedVect1[dataNum], ionBackFlowCorrectedErrorVect1[dataNum], ionBackFlowCorrectedVect1_old[dataNum], ionBackFlowCorrectedErrorVect1_old[dataNum];
	LoadIbfData(gasName, dataNum, hvMeshListIBF1, hvDriftListIBF1, ionBackFlowCorrectedVect1, ionBackFlowCorrectedErrorVect1, ionBackFlowCorrectedVect1_old, ionBackFlowCorrectedErrorVect1_old);
	
	
	Double_t ratioListData[dataNum];
	for (int i = 0; i < dataNum; i++) {ratioListData[i] = hvMeshListIBF1[i]/(hvDriftListIBF1[i]-hvMeshListIBF1[i]) * (ddrift-damp)/damp;}
	
	TGraphErrors* grData1 = new TGraphErrors(dataNum, hvMeshListIBF1, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
	grData1->SetMarkerStyle(20);
	grData1->SetMarkerSize(0.3);
	grData1->SetMarkerColor(4);
	grData1->Draw("LP same");
	TGraphErrors* grData2 = new TGraphErrors(dataNum, hvMeshListIBF1, ionBackFlowCorrectedVect1_old, 0, ionBackFlowCorrectedErrorVect1_old );
	grData2->SetMarkerStyle(20);
	grData2->SetMarkerSize(0.3);
	grData2->SetMarkerColor(5);
	grData2->Draw("LP same");
	
	TLegend* legend = new TLegend(0.5,0.65,0.9,0.9);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legend->AddEntry(grData1,"Data","lp");
	legend->AddEntry(grData2,"Data old","lp");
	legend->AddEntry(grSim1,"Simulation (counting IBF)","lp");
	legend->AddEntry(grSim2, "Simulation (induced charge)", "lp");
	legend->AddEntry(grSim3, "Simulation (ratio of integrated currents)", "lp");
	legend->Draw();
	
	
	c5->cd(2);
	TGraphErrors* grSimRatio1 = new TGraphErrors(num, hvRatioList, ibfList, 0, ibfErrorList);
	grSimRatio1->SetTitle("IBF = f(E ratio)");
	grSimRatio1->GetXaxis()->SetTitle( "E_{amp}/E_{drift}" );
	grSimRatio1->GetYaxis()->SetTitle( "IBF (%)" );
	//grSimRatio1->GetXaxis()->SetLimits(0, 4);  // along X axis
	grSimRatio1->GetHistogram()->SetMinimum(0.);   // along Y axis
	grSimRatio1->GetHistogram()->SetMaximum(4.);   // along Y axis
	grSimRatio1->SetMarkerStyle(20);
	grSimRatio1->SetMarkerSize(0.3);
	grSimRatio1->SetMarkerColor(1);
	//grSimRatio1->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num]+5);
	grSimRatio1->Draw("ALP");
	
	TGraphErrors* grSimRatio2 = new TGraphErrors(num, hvRatioList, ibfList2, 0, ibfErrorList);
	grSimRatio2->SetTitle("IBF = f(E ratio)");
	grSimRatio2->GetXaxis()->SetTitle( "E_{amp}/E_{drift}" );
	grSimRatio2->GetYaxis()->SetTitle( "IBF (%)" );
	//grSimRatio2->GetXaxis()->SetLimits(0, 4);  // along X axis
	grSimRatio2->GetHistogram()->SetMinimum(0.);   // along Y axis
	grSimRatio2->GetHistogram()->SetMaximum(4.);   // along Y axis
	grSimRatio2->SetMarkerStyle(20);
	grSimRatio2->SetMarkerSize(0.3);
	grSimRatio2->SetMarkerColor(2);
	//grSimRatio2->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num]+5);
	grSimRatio2->Draw("LP same");
	
	TGraph* grSimRatio3 = new TGraph(num, hvRatioList, ibfList3);
	grSimRatio3->SetTitle("IBF = f(E ratio)");
	grSimRatio3->GetXaxis()->SetTitle( "E_{amp}/E_{drift}" );
	grSimRatio3->GetYaxis()->SetTitle( "IBF (%)" );
	//grSimRatio3->GetXaxis()->SetLimits(0, 4);  // along X axis
	grSimRatio3->GetHistogram()->SetMinimum(0.);   // along Y axis
	grSimRatio3->GetHistogram()->SetMaximum(4.);   // along Y axis
	grSimRatio3->SetMarkerStyle(20);
	grSimRatio3->SetMarkerSize(0.3);
	grSimRatio3->SetMarkerColor(3);
	//grSimRatio2->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num]+5);
	grSimRatio3->Draw("LP same");
	
	
	TGraphErrors* grDataRatio1 = new TGraphErrors(dataNum, ratioListData, ionBackFlowCorrectedVect1, 0, ionBackFlowCorrectedErrorVect1 );
	grDataRatio1->SetMarkerStyle(20);
	grDataRatio1->SetMarkerSize(0.3);
	grDataRatio1->SetMarkerColor(4);
	grDataRatio1->Draw("LP same");
	TGraphErrors* grDataRatio2 = new TGraphErrors(dataNum, ratioListData, ionBackFlowCorrectedVect1_old, 0, ionBackFlowCorrectedErrorVect1_old );
	grDataRatio2->SetMarkerStyle(20);
	grDataRatio2->SetMarkerSize(0.3);
	grDataRatio2->SetMarkerColor(5);
	grDataRatio2->Draw("LP same");
	
	TLegend* legendRatio = new TLegend(0.5,0.65,0.9,0.9);
	//legendRatio->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legendRatio->AddEntry(grDataRatio1, "Data", "lp");
	legendRatio->AddEntry(grDataRatio2, "Data old", "lp");
	legendRatio->AddEntry(grSimRatio1,"Simulation (counting IBF)", "lp");
	legendRatio->AddEntry(grSimRatio2, "Simulation (induced charge)", "lp");
	legendRatio->AddEntry(grSimRatio3, "Simulation (ratio of integrated currents)", "lp");
	legendRatio->Draw();
	
	c5->SaveAs(Form("Figures/IBFCurve-%s-model%d.pdf", gasName.c_str(), modelNum));
	
	
	time_t t1 = time(NULL);
	//PrintTime(t0, t1);
	
	return 0;
}


