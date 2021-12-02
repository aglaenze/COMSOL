#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <TColor.h>
#include <TCanvas.h>

#include "Include/Utils.C"

using namespace std;
/*
 I simulate root files with 100 events, and if I want to have them as one file, I have to add current
 ! but there's an overlapping period (ionDelay) where currents need to be added
 */

Double_t Square(Double_t x) {
    return x*x;
}

Double_t FitFunction(Double_t* x, Double_t* par ) { //(x, alpha, n sigma, mu)
    return 1-exp(-x[0]*par[0]*par[1]);
}

int PhotonConversionRatio() {
    
    //______________________
    // variables
    //string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //string gasName = "air";
    //____________________
    
    //double hAbove = 10;      //cm
    double hAbove = 0;      //cm
    double lDetector = 12;   // size of active area in x, z (in cm)
    double hDetector = 3;    // height of the detector in z (in cm)
    
    //double zMin = 2.265;    // distance from the top of the detector to the closest electrode of interest
    //double zMax = 2.987;    // distance from the top of the detector to the furthest electrode of interest
    
    /*
    double zMin = 1.143;
    double zMax = 2.193;
     */
    double zMin = 2.487;
    double zMax = 2.987;
    
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.07);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetTitleSize(.08);
    gStyle->SetTextSize(.08);
    gStyle->SetLabelSize(.05, "XY");
    
    const char* filename = Form("rootFiles/%s/ConversionNumber.root", gasName.c_str());
    TFile* f = new TFile(filename, "READ");
    
    TH1* hConversion = (TH1*)f->Get("hConversion");
    
    string gasTitle = gasName;
    if (gasName == "Ar-iC4H10") gasTitle = "Ar-isobutane (95/5)";
    else if (gasName == "air") gasTitle = "common air";
    TCanvas* cv = new TCanvas("conv", "Proportion of photons that converted to electrons", 300, 200);
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.1);
    hConversion->SetXTitle("length traversed (cm)");
    hConversion->SetYTitle("P_{#gamma}");
    hConversion->SetTitle(Form("Proportion of photons that converted to electrons in %s", gasTitle.c_str()));
    hConversion->SetMinimum(0);
    hConversion->Draw("hist");
    
    
    // Now find the theoretical ratio, using the cross section of Ar and gamma at 5.9keV
    double ddrift = 0.;
    double ddriftMin = 0.3;
    double zStep = 0.2;
    int nStep = 50;
    double ddriftMax = ddriftMin+nStep*zStep; // cm
    
    // properties of Argon
    double densityMass = 1.784e-3;      // g/cm3
    double densityMol = 39.95;          // g/mol
    double AvNumber = 6.0221409e+23;    // mol-1
    
    double density = 0.95*densityMass/densityMol*AvNumber;  // factor 0.95 for the proportion of Argon in the mixture
    
    cout << "Argon density in the mixture = " << density << " cm-3" << endl;
    
    double crossSection = 1.6e-20;  // density x cross section has to be in cm-1    // A CHANGER !! TROUVER CETTE VALEUR POUR 5.9 keV
    //For gases the one-photon absorption cross-section σ1 is typically of the order of 10−17cm2[20], whereas the two-photon and the three-photon cross-sections are of the order of σ2 = W/F 2 ∼ 10−50cm4s and σ3 = W/F 3 ∼ 10−83cm6s2, respectively.
    // PDG chap 33.4.5, not an exact value but at least it's compatible
    
    /*
     TH1F* hTh = new TH1F("hTh", "Proportion of photons that converted (th)", nStep, ddriftMin-zStep/2, ddriftMin+nStep*zStep-zStep/2);
     for (int k = 0; k<nStep; k++) {
     ddrift = ddriftMin+k*zStep; // cm
     hTh->Fill(ddrift, 1-exp(-density*crossSection*ddrift));
     }
     hTh->SetLineColor(2);
     hTh->Draw("hist same");
     */
    TF1* fTh = new TF1("FitFunction", FitFunction, ddriftMin, ddriftMax, 2);
    fTh->FixParameter(0, density);
    fTh->SetParameter(1, crossSection);
    hConversion->Fit(fTh);
    fTh->Draw("same");
    
    double xText = ddriftMin + (ddriftMax-ddriftMin)*0.4;
    TLatex* txt = new TLatex(xText,0.3,Form("#bf{#sigma_{#gamma-Ar} = %.3g cm^{-2}}", fTh->GetParameter(1)));
    txt->Draw();
    
    cv->SaveAs(Form("Figures/PhotonConversionRatio-%s.pdf", gasName.c_str()));
    
    /*
     // Compute the proportion of photons that converted within a gap
     double z0 = 1.5;          // cm, correspond à l'espace entre le toit du détecteur (position de la source) et la dérive
     double zMin = z0 + 1.122;    // cm
     double zMax = z0 + 1.844;    // cm
     
     int iBinMin = hConversion->FindBin(zMin);
     int iBinMax = hConversion->FindBin(zMax);
     
     cout << endl << "Proportion of photons that converted between zmin = " << zMin << " and zmax = " << zMax << " is " << hConversion->GetBinContent(iBinMax) - hConversion->GetBinContent(iBinMin) << endl;
     */
    
    // Compute the fraction of photons that converted in the space of interest
    
    // 1) Get the fraction of photons that did not convert in the air above the detector
    const char* filenameAir = "rootFiles/air/ConversionNumber.root";
    TFile* fAir = new TFile(filenameAir, "READ");
    TH1* hConversionAir = (TH1*)fAir->Get("hConversion");
    
    double fraction = 0;
    double thetaMax = TMath::ATan(lDetector/2 /(hDetector+hAbove));
    //double thetaMax = 0.25*TMath::Pi();
    cout << "thetaMax = " << thetaMax << endl;
    int nStepTheta = 100;
    int nStepZ = 100;
    double stepSizeZ = (zMax-zMin)/nStepZ;
    double stepSizeTheta = thetaMax/nStepTheta;
    
    // Loop starts here over theta
    for (int i = 0; i < nStepTheta; i++) {
        double theta = i*stepSizeTheta;
        int iBinAir = hConversionAir->FindBin(hAbove/TMath::Cos(theta));
        double propAir = 1-hConversionAir->GetBinContent(iBinAir);   // proba that photons did not convert in the air
        if (iBinAir > hConversionAir->GetNbinsX()-1) propAir = 0;
        //cout << "propAir = " << propAir << endl;
        // Loop starts over z
        for (int j = 0; j < nStepZ; j++) {
            double z = zMin + (j+1)*stepSizeZ;
            double r = z/TMath::Cos(theta);
            double zBefore = zMin + j*stepSizeZ;
            double rBefore = zBefore/TMath::Cos(theta);
            int iBin = hConversion->FindBin(r);
            int iBinBefore = hConversion->FindBin(rBefore);
            double propConv = hConversion->GetBinContent(iBin)-hConversion->GetBinContent(iBinBefore);
            if (iBin > hConversion->GetNbinsX()-1) propConv = 0;
            double probaFin = propAir*propConv;

            fraction += propConv * stepSizeTheta/TMath::Pi();
            //cout << "fraction = " << fraction << endl;
        }   // end of loop on z
    }   // end of loop on theta
    
    cout << endl << "Proportion of photons that converted between zmin = " << zMin << " and zmax = " << zMax << " and between thetaMin = 0 and thetaMax = " << thetaMax/TMath::Pi() << " pi is " << fraction*100 << "%" << endl;
    

    return 0;
}


