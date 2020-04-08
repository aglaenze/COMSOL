#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "_Utils.C"

/*
 I simulate root files with 100 events, and if I want to have them as one file, I have to add current
 ! but there's an overlapping period (ionDelay) where currents need to be added
 */

int AddIntegralBranch() {
    
    int modelNum = 1;
    TString path = Form("rootFiles/Ar-iC4H10/model%d/", modelNum);
    TString fileName = path + "signal-360-560-1.root";
    TFile *f = new TFile(fileName,"update");

    int electrodeNum = GetElectrodeNum(modelNum);
    std::cout << "there are " << electrodeNum << " electrodes" << std::endl;
    for (int k = 0; k<electrodeNum; k++) {
        TTree *T = (TTree*)f->Get(Form("tSignal_%d", k+2));
        int nEntries = T->GetEntries();
        double fctInt = 0, fceInt = 0, fciInt = 0;
        TBranch *bStInt = T->Branch("totalCurrentInt",&fctInt,"totalCurrentInt/D");
        TBranch *bSeInt = T->Branch("electronCurrentInt",&fceInt,"electronCurrentInt/D");
        TBranch *bSiInt = T->Branch("ionCurrentInt",&fciInt,"ion/D");
        double fct, fce, fci;
        T->SetBranchAddress("totalCurrent",&fct);
        T->SetBranchAddress("electronCurrent",&fce);
        T->SetBranchAddress("ionCurrent",&fci);
        for (int l = 0; l<nEntries; l++) {
            T->GetEntry(l);
            fctInt+=fct;
            fceInt+=fce;
            fciInt+=fci;
            bStInt->Fill();
            bSeInt->Fill();
            bSiInt->Fill();
        }
        T->Write();
    }
    return 0;

}


