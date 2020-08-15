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
