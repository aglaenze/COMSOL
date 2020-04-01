#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

/*
 I simulate root files with 20 events, and if I want to have them as one file, I have to add current
 ! but there's an overlapping period (ionDelay) where currents need to be added
 */

int AddSignalTrees() {
    
    TString path = "rootFiles/Ar-iC4H10/model1/";

    TString outputName = path + "feSignal-340-540.root";
    TFile* ff = new TFile(outputName, "RECREATE");
    
    const int electrodeNum = 3;
    const int nEvents = 20;
    const double tStep = 0.1;   //ns
    const double rate = 6.e7;               // number of events per s (note that t units are ns here, we'll need a conversion factor)
    const double timespace = 1./rate*1.e9;    // in ns
    double signalDuration = nEvents * timespace;
    
    const int numberOfFiles = 2;
    
    for (int k = 0; k < electrodeNum; k++) {
        // new TTree
        TTree *tSignalNew = new TTree(Form("tSignal_%d",k+2),"Currents");
        Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
        tSignalNew->Branch("time", &ft, "time/D");
        tSignalNew->Branch("totalCurrent", &fct, "totalCurrent/D");
        tSignalNew->Branch("electronCurrent", &fce, "electronCurrent/D");
        tSignalNew->Branch("ionCurrent", &fci, "ionCurrent/D");
        
        int nLimit = 0;
        
        for (int i = 1; i<numberOfFiles+1; i++) {
            TString fileName = path + Form("feSignal-340-540-%d.root", i);
            std::cout << "fileName = " << fileName << std::endl;
            TFile* f = TFile::Open(fileName, "READ");
            
            /*
            TString fileName2;
            TFile* f2;
            if (i<numberOfFiles) {
                fileName2 = path + Form("feSignal-340-540-%d.root", i+1);
                std::cout << "fileName2 = " << fileName2 << std::endl;
                f2 = TFile::Open(fileName2, "READ");
            }
             */
        
            // data from first tree
            Double_t tt = 0., ct = 0., ce = 0., ci = 0.;
            TTree* tSignal = (TTree*)f->Get(Form("tSignal_%d", k+2));
            tSignal->SetBranchAddress("time", &tt);
            tSignal->SetBranchAddress("totalCurrent", &ct);
            tSignal->SetBranchAddress("electronCurrent", &ce);
            tSignal->SetBranchAddress("ionCurrent", &ci);
            
            /*
            // data from next overlaping tree
            Double_t tt2 = 0., ct2 = 0., ce2 = 0., ci2 = 0.;
            TTree* tSignal2;
            if (i<numberOfFiles) {
                tSignal2 = (TTree*)f2->Get(Form("tSignal_%d", k+2));
                tSignal2->SetBranchAddress("time", &tt2);
                tSignal2->SetBranchAddress("totalCurrent", &ct2);
                tSignal2->SetBranchAddress("electronCurrent", &ce2);
                tSignal2->SetBranchAddress("ionCurrent", &ci2);
            }
             */

            int n = tSignal->GetEntries();
            
            for (int l = nLimit; l < n; l++) {
                tSignal->GetEntry(l);
                ft = (i-1)*signalDuration + tt;
                fct = ct;
                fce = ce;
                fci = ci;
                //return 0;
                if (ft == i*signalDuration) {
                    nLimit = l;
                    std::cout << "nLimit = " << nLimit << std::endl;}
                /*
                if (ft > i*signalDuration && i<numberOfFiles) {
                    tSignal2->GetEntry(l-nLimit);
                    fct += ct2;
                    fce += ce2;
                    fci += ci2;
                }
                 */
                //return 0;
                tSignalNew->Fill();
                //return 0;
            }
            f->Close();
            //f2->Close();
            //return 0;
        }
        //return 0;
        tSignalNew->Write();
    }
    //return 0;
    ff->Close();

    return 0;

}


