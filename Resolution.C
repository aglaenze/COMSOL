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
	
	vector<Double_t> ibf = {1, 3};
	vector<Double_t> ibfError = {0.5, 0.5};
	
	vector<Double_t> resolution = {10, 15};
	vector<Double_t> resolutionError = {2, 2};
	

	TCanvas* cv = new TCanvas("cv", "cv", 800, 800);
	TGraphErrors* tResolution = CreateTGraph( ibf, resolution, ibfError, resolutionError );
	tResolution->SetTitle("#sigma(E)/E = f(IBF)");
	tResolution->GetYaxis()->SetTitle("#sigma(E)/E (%)");
	tResolution->GetXaxis()->SetTitle("IBF (%)");
	tResolution->Draw();
	cv->SaveAs(Form("Figures/model%d/Resolution-IBF.pdf", modelNum));
	
	return 0;
}
