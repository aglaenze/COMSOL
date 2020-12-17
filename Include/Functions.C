#include <iostream>
#include <fstream>
#include <sstream>
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


Double_t Square(Double_t x ) {return x*x;}

Double_t FitGauss(Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*TMath::Gaus( x[0], par[0], par[1]); }

//____________________________________________
Double_t FitFunctionExp(Double_t* x, Double_t* par) {
	return TMath::Exp( par[0] + par[1]*x[0] ); }

Double_t FitLandau(Double_t* x, Double_t* par) { //Double_t TMath::Landau	(Double_t x, Double_t mu = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
	return  par[2]*TMath::Landau( x[0], par[0], par[1]); }

Double_t CrystalBall(Double_t x, Double_t mean, Double_t sigma, Double_t alpha, Double_t n)
{
	
	Double_t t = (x-mean)/sigma;
	if( alpha < 0 ) t *= -1.0;
	
	alpha = fabs( alpha );
	if( t >= -alpha ) return TMath::Exp( -Square( t )/2 );
	else {
		
		Double_t a = TMath::Power( n/alpha, n )*TMath::Exp( -Square( alpha )/2 );
		Double_t b = n/alpha - alpha;
		return a/TMath::Power( b - t, n );
	}
}

Double_t FitFunctionCrystalBall( Double_t* x, Double_t* par ) {
	return par[0]*CrystalBall( x[0], par[1], par[2], par[3], par[4] ); }


TF1* GetFitCurve(TH1F* h, bool gauss = true) {
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	std::cout << "xMax = " << xMax << std::endl;
	std::cout << "maximum = " << h->GetMaximum() << std::endl;
	
	Int_t fitRangeMin = 0;
	Int_t fitRangeMax = xMax + 0.5 * h->GetRMS();
	if (gauss) {
		fitRangeMin = xMax - 1.1 * h->GetRMS();
		fitRangeMax = xMax + 1.5*h->GetRMS();
	}
	if (fitRangeMin < 0) fitRangeMin = 0;
	
	TF1* f;
	if (gauss) f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
	else f = new TF1( "FitFunction", FitLandau, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("Mean", "Sigma", "Amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	
	//std::cout << "\n\nh->GetRMS() = " << h->GetRMS() << std::endl;
	//std::cout << "\n\nh->GetMaximum() = " << h->GetMaximum() << std::endl;
	
	if (!gauss) {
		f->SetParLimits(0, 0.5*xMax, 2*xMax);
		//f->SetParLimits(1,0, 0.1*h->GetRMS());
		//f->FixParameter(2, 1);
		//f->SetParLimits(2, 0.5*h->GetMaximum(), 5*h->GetMaximum());
		//f->SetParameter(2,4*h->GetMaximum());
	}
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}

TF1* GetFitIbf(TH1F* h, bool gauss = true) {
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	std::cout << "xMax = " << xMax << std::endl;
	std::cout << "maximum = " << h->GetMaximum() << std::endl;
	
	Int_t fitRangeMin = xMax - h->GetRMS();
	//Int_t fitRangeMin = 0;
	Int_t fitRangeMax = xMax + 1*h->GetRMS();
	
	TF1* f;
	if (gauss) f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
	else f = new TF1( "FitFunction", FitLandau, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("Mean", "Sigma", "Amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	f->SetParLimits(0,0,10);
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}

