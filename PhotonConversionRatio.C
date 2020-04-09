#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <TColor.h>
#include <TCanvas.h>

#include "_Utils.C"

/*
 I simulate root files with 100 events, and if I want to have them as one file, I have to add current
 ! but there's an overlapping period (ionDelay) where currents need to be added
 */

int PhotonConversionRatio() {
    
    //______________________
    // variables
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //____________________
    
    
    const char* filename = Form("rootFiles/%s/ConversionNumber.root", gasName.c_str());
    TFile* f = new TFile(filename, "READ");
    
    TH1* hConversion = (TH1*)f->Get("hConversion");
    
    gStyle->SetOptStat(0);
    TCanvas* cv = new TCanvas("conv", "Proportion of photons that converted to electrons", 600, 300);
    hConversion->SetXTitle("length traversed (cm)");
    hConversion->Draw("hist");
    
    
    // Now find the theoretical ratio, using the cross section of Ar and gamma at 5.9keV
    double ddrift = 0.;
    double ddriftMin = 0.3;
    double zStep = 0.2;
    int nStep = 50;
    
    TH1F* hTh = new TH1F("hTh", "Proportion of photons that converted (th)", nStep, ddriftMin-zStep/2, ddriftMin+nStep*zStep-zStep/2);
    
    // properties of Argon
    double densityMass = 1.784e-3;      // g/cm3
    double densityMol = 39.95;          // g/mol
    double AvNumber = 6.0221409e+23;    // mol-1

    double density = 0.95*densityMass/densityMol*AvNumber;  // factor 0.95 for the proportion of Argon in the mixture
    
    std::cout << "Argon density in the mixture = " << density << " cm-3" << std::endl;
    
    double crossSection = 1.6e-20;  // density x cross section has to be in cm-1    // A CHANGER !! TROUVER CETTE VALEUR POUR 5.9 keV
    
    for (int k = 0; k<nStep; k++) {
        ddrift = ddriftMin+k*zStep; // cm
        hTh->Fill(ddrift, 1-exp(-density*crossSection*ddrift));
    }
    
    hTh->SetLineColor(2);
    hTh->Draw("hist same");
    
    
    cv->SaveAs("Figures/PhotonConversionRatio.pdf");
    return 0;
}


