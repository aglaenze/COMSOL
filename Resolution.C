#include <iostream>
#include <fstream>
#include <cmath>
#include <dirent.h>
#include <string.h>
#include <stdio.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TString.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

using namespace std;

//____________________________________________
TGraphErrors* CreateTGraph( Int_t size, const Double_t* x, const Double_t* y, const Double_t* xErr, const Double_t* yErr )
{
	TGraphErrors* tg = new TGraphErrors();
	for( Int_t i = 0; i < size; ++i )
	{
		tg->SetPoint( i, x[i], y[i] );
		tg->SetPointError( i, xErr[i], yErr[i] );
	}
	
	return tg;
}

//____________________________________________
TGraphErrors* CreateTGraph( const std::vector<Double_t>& x, const std::vector<Double_t>& y, const std::vector<Double_t>& xErr, const std::vector<Double_t>& yErr )
{ return CreateTGraph( x.size(), &x[0], &y[0], &xErr[0], &yErr[0] ); }

int Resolution() {
	//std::map <std::string, int, cicompare> electrodeMap;
	
	int modelNum = 10;
	
	vector<Double_t> ibf = {0.0156601, 0.0025, 0.0360657, 0.0175};
	vector<Double_t> ibfError = {0.0189029, 0.0051282, 0.0130104, 0.0520183};
	
	vector<Double_t> resolution = {0.465261, 0.513275, 1.13225, 0.531707};
	vector<Double_t> resolutionError = {0.0758699, 0.062831, 0.0130104, 0.0520183};
	
	for (int k = 0; k<resolution.size(); k++) {resolution[k] = resolution[k]*100; resolutionError[k] = resolutionError[k]*100;}

	TCanvas* cv = new TCanvas("cv", "cv", 800, 800);
	TGraphErrors* tResolution = CreateTGraph( ibf, resolution, ibfError, resolutionError );
	tResolution->SetTitle("#sigma(E)/E = f(IBF)");
	tResolution->GetYaxis()->SetTitle("#sigma(E)/E (%)");
	tResolution->GetXaxis()->SetTitle("IBF (%)");
	tResolution->Draw();
	cv->SaveAs(Form("Figures/model%d/Resolution-IBF.pdf", modelNum));
	
	return 0;
}

/*
330-410-530-680-800.pdf has been
Mean Gain = 666.125 pm 59.9042
Sigma Gain = 309.922 pm 42.159
Resolution = 0.465261 pm 0.0758699
IBF = 0.0156601 pm 0.0189029

330 410 530 730 850
Mean Gain = 779.595 pm 53.024
Sigma Gain = 400.147 pm 40.726
Resolution = 0.513275 pm 0.062831
IBF = 0.0025 pm 0.0051282

330 410 530 830 950
Mean Gain = 1801.77 pm 562.447
Sigma Gain = 2040.06 pm 239.132
Resolution = 1.13225 pm 0.377544
IBF = 0.0360657 pm 0.0130104

350 430 530 830 950
Mean Gain = 4107.45 pm 254.11
Sigma Gain = 2183.96 pm 197.33
Resolution = 0.531707 pm 0.0582244
IBF = 0.0175 pm 0.0520183
*/
