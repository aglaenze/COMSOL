#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>


int RepareFile() {
    
    TString path = "rootFiles/Ar-CO2/model4/";
    //TString fileName = path + "ibf_420V.root";
    TString fileName = path + "gain_350V.root";

// Open broken file to get the histogram
    TFile* f = TFile::Open(fileName, "READ");
    //TH1F* h = (TH1F*)f->Get("hFeElectrons");
    //TH1F* h = (TH1F*)f->Get("hIbf");
    TH1F* h = (TH1F*)f->Get("hElectrons");

// Open a new file to copy the histogram there and close it
    TFile* ff = new TFile(fileName, "RECREATE");
    h->Write();
    ff->Close();
    
    f->Close();
    //f->Delete();

    return 0;

}


