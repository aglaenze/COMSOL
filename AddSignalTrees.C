#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

/*
 I simulate root files with 100 events, and if I want to have them as one file, I have to add current
 ! but there's an overlapping period (ionDelay) where currents need to be added
 */

int AddSignalTrees() {
    
    TString path = "rootFiles/Ar-iC4H10/model1/";
    TString outputName = path + "feSignal-340-540.root";
    TFile* fOut = new TFile(outputName, "RECREATE");
    
    const int numberOfFiles = 12;
    
    TTree *tGain = new TTree("tGain","Gain");
    TTree *tZstart = new TTree("tZstart","Departure z of electrons");
    Int_t gain = 0;
    Double_t zStart = 0.;
    tGain->Branch("gain", &gain, "gain/I");
    tZstart->Branch("zStart", &zStart, "zStart/D");
    
    for (int i = 1; i<numberOfFiles+1; i++) {
        TString fileName = path + Form("feSignal-340-540-%d.root", i);
        TFile* fIn = TFile::Open(fileName, "READ");
        //Int_t ne2 = 0;
        //Double_t ze1 = 0;
        TTree* tNe = (TTree*)fIn->Get("tGain");
        //tNe->SetBranchAddress("electronNumber", &ne2);
        tNe->SetBranchAddress("electronNumber", &gain);
        
        int nNe = tNe->GetEntries();
        for (int l = 0; l<nNe; l++) {
            tNe->GetEntry(l);
            //gain = ne2;
            tGain->Fill();
        }
       
        TTree* tZe = (TTree*)fIn->Get("tZstart");
        //tZe->SetBranchAddress("zStart", &ze1);
        tZe->SetBranchAddress("zStart", &zStart);
        
        int nZe = tZe->GetEntries();
        for (int l = 0; l<nZe; l++) {
            tZe->GetEntry(l);
            tZstart->Fill();
        }
         
        fIn->Close();
    }
    fOut->cd();
    tGain->Write();
    tZstart->Write();
    fOut->Close();
    

    return 0;
    
    //__________________//
    // à partir de là ça ne marche pas
    
    const int electrodeNum = 3;
    const int nEvents = 100;
    const double tStep = 0.1;   //ns
    const double rate = 6.e7;               // number of events per s (note that t units are ns here, we'll need a conversion factor)
    const double timespace = 1./rate*1.e9;    // in ns
    double signalDuration = nEvents * timespace;
    
    std::vector<TTree*> treeVector;
    
    
    for (int k = 0; k < electrodeNum; k++) {
        // new TTree
        TTree *tSignalNew = new TTree(Form("tSignal_%d",k+2),"Currents");
        Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
        tSignalNew->Branch("time", &ft, "time/D");
        tSignalNew->Branch("totalCurrent", &fct, "totalCurrent/D");
        tSignalNew->Branch("electronCurrent", &fce, "electronCurrent/D");
        tSignalNew->Branch("ionCurrent", &fci, "ionCurrent/D");
        
        int nLimit = 0;
        double signalDurationT = 0;
        double tPast = 0.;
        double tRef = 0;
        
        for (int i = 1; i<numberOfFiles+1; i++) {
            signalDurationT = int(10*signalDuration*i);
            signalDurationT /= 10.;
            std::cout << "i*signalDuration = " << signalDurationT << std::endl;
            
            TString fileName = path + Form("feSignal-340-540-%d.root", i);
            std::cout << "fileName = " << fileName << std::endl;
            TFile* f = TFile::Open(fileName, "READ");

            TString fileName2;
            TFile* f2;
            if (i<numberOfFiles) {
                fileName2 = path + Form("feSignal-340-540-%d.root", i+1);
                std::cout << "fileName2 = " << fileName2 << std::endl;
                f2 = TFile::Open(fileName2, "READ");
            }
        
            // data from first tree
            Double_t tt = 0., ct = 0., ce = 0., ci = 0.;

            TTree* tSignal = (TTree*)f->Get(Form("tSignal_%d", k+2));
            tSignal->SetBranchAddress("time", &tt);
            tSignal->SetBranchAddress("totalCurrent", &ct);
            tSignal->SetBranchAddress("electronCurrent", &ce);
            tSignal->SetBranchAddress("ionCurrent", &ci);

            
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
            
            //int n = 1000;
            int n = tSignal->GetEntries();
            std::cout << "tPast = " << tPast << std::endl;
            
            //Double_t tt = 0., ct = 0., ce = 0., ci = 0.;
            for (int l = nLimit; l < n; l++) {
                //if (l-nLimit < 0) return 2;
                tSignal->GetEntry(l);
                ft = tPast + tt;
                //if (ft-tRef != tStep) {std::cout << "problem in time, dt = " << ft-tRef << std::endl;}
                //std::cout << "ft = " << ft << std::endl;
                tRef = tPast + tt;
                fct = ct;
                fce = ce;
                fci = ci;
                //return 0;
        
                if (ft == signalDurationT) {
                    nLimit = l;
                    std::cout << "\n\nnLimit = " << nLimit << std::endl;
                }
                 
                
                if (ft > signalDurationT && i<numberOfFiles) {
                    tSignal2->GetEntry(l-nLimit);
                    fct += ct2;
                    fce += ce2;
                    fci += ci2;
                }

                tSignalNew->Fill();
            }   // end for l < nEntries
            f->Close();
            f2->Close();
            tPast+=signalDurationT;
                //return 0;
        }   // end for i<numberOfFiles
        //return 0;
        treeVector.push_back(tSignalNew);
        //tSignalNew->Write();
        
    } // end for k < electrodeNum
    
    for (int i = 0; i<treeVector.size(); i++) treeVector[i]->Write();
    fOut->Close();
    
    return 0;

}


