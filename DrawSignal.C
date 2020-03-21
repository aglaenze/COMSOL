#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "_Utils.C"
#include "_Data.C"



int DrawSignal() {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    //____________________

    time_t t0 = time(NULL);
    gStyle->SetOptStat(0);
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    TFile fFe(Form("rootFiles/%s/spectrum_Fe55.root", gasName.c_str()));
    TH1F* hFe = (TH1F*)fFe.Get("hElectrons");
    Int_t nPrimaryTh = GetPrimary(gasName);
     
    // Get number of files to look at
    Int_t num = 1;
    //Int_t num = GetNumberOfFiles(path, "signal");
    /*Int_t num2 = int(num/2.);
    if (num/2.> num2) num2+=1;
    TCanvas* c2 = new TCanvas("c2");
    c2->Divide(2, num2);
    */
    Int_t electrodeNum = 0;
    if (modelNum == 1) electrodeNum = 3;
    else if (modelNum >= 2 && modelNum < 5) electrodeNum = 4;
    else if (modelNum >= 5 && modelNum < 10) electrodeNum = 5;
    
    TCanvas* cvf = new TCanvas("cvf");
    cvf->Divide(electrodeNum, num);
    
    Double_t hvMeshList[num];
    Double_t gainList[num];
    Double_t gainErrorList[num];
    
    /*
    for (unsigned int k = 0; k < num; ++k) {
        Int_t hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvGemDown = 0, hvGemUp = 0, hvDrift = 0;
        TString signalFileName, fOutputName;
        if (modelNum == 1) {
            hvMesh = 340+20*k;
            hvDrift = 540+20*k;
            signalFileName = path + Form("signal-%d-%d.root", hvMesh, hvDrift);
            fOutputName = path + Form("amp-%d-%d.root", hvMesh, hvDrift);
            //
            signalFileName = path + "signal-380-580.root";
            fOutputName = path + "amp-380-580.root";
        }
        else if (modelNum >= 2 && modelNum < 5) {
            hvDmDown = 300;
            hvDmUp = 600;
            hvDrift = 800;
            signalFileName = path + Form("signal-%d-%d-%d.root", hvDmDown, hvDmUp, hvDrift);
            fOutputName = path + Form("amp-%d-%d-%d.root", hvDmDown, hvDmUp, hvDrift);
        }
        else if (modelNum >= 5 && modelNum < 8) {
            hvMesh = 300;
            hvDmDown = 400;
            hvDmUp = 700;
            hvDrift = 900;
            signalFileName = path + Form("signal-%d-%d-%d-%d.root", hvMesh, hvDmDown, hvDmUp, hvDrift);
            fOutputName = path + Form("amp-%d-%d-%d-%d.root", hvMesh, hvDmDown, hvDmUp, hvDrift);
        }
        else if (modelNum >= 8 && modelNum < 14) {
            hvMesh = 350;
            hvDmDown = 430;
            hvDmUp = 630;
            hvDrift = 750;
            signalFileName = path + Form("signal-%d-%d-%d-%d.root", hvMesh, hvGemDown, hvGemUp, hvDrift);
            fOutputName = path + Form("amp-%d-%d-%d-%d.root", hvMesh, hvGemDown, hvGemUp, hvDrift);
        }
     */
        // import data
        //TFile* fSignal = TFile::Open(signalFileName, "READ");
        //TFile* fAmp = new TFile(fOutputName, "RECREATE");
    
    TFile* fSignal = TFile::Open("rootFiles/Ar-iC4H10/model1/signal-340-540.root", "READ");
    TFile* fAmp = new TFile("rootFiles/Ar-iC4H10/model1/charges-340-540.root", "RECREATE");
    
    std::map <std::string, int> electrode;
    if (modelNum==1) {
        electrode["mesh"] = 2;
        electrode["drift"] = 3;
        electrode["pad"] = 4;
    }

  
        for (int j = 0; j < electrodeNum; j++) {
        //for (int j = 0; j < 1; j++) {
            
            /*
            TCanvas* c2 = new TCanvas("c2", "c2", 1500, 500);
            c2->Divide(3);
            
            TH1F* hIon = new TH1F(Form("hIon_V%d", j+2), Form("hIon_V%d", j+2), 100, 0, 0.2);
            TH1F* hElectron = new TH1F(Form("hElectron_V%d", j+2), Form("hElectron_V%d", j+2), 300, 0, 3000);
            TH1F* hTotal = new TH1F(Form("hTotal_V%d", j+2), Form("hTotal_V%d", j+2), 300, 0, 3000);
             */
            
            TTree *tCharges = new TTree(Form("tCharges_V%d", j+2),Form("tCharges_V%d", j+2));
            Double_t nIons = 0, nElectrons = 0, nTotal = 0;
            tCharges->Branch("nIons", &nIons, "nIons/D");
            tCharges->Branch("nElectrons", &nElectrons, "nElectrons/D");
            tCharges->Branch("nTotal", &nTotal, "nTotal/D");
   
            
            Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
            TTree* tSignal = (TTree*)fSignal->Get(Form("tSignal_V%d", j+2));
            tSignal->SetBranchAddress("time", &ft);
            tSignal->SetBranchAddress("totalCurrent", &fct);
            tSignal->SetBranchAddress("electronCurrent", &fce);
            tSignal->SetBranchAddress("ionCurrent", &fci);
            
            //Double_t chargeNumberIon = 0., chargeNumberElectron = 0., chargeNumberTotal = 0.;
            double thresholdPos = tSignal->GetMaximum("ionCurrent");
            double thresholdNeg = tSignal->GetMinimum("ionCurrent");
            double threshold = 0;
            if ( abs(thresholdNeg) > abs(thresholdPos)) threshold = abs(thresholdNeg);
            else threshold = abs(thresholdPos);
            threshold = 0.2;
            std::cout << "threshold = " << threshold << std::endl;
            
            
            int nEvents = 0;
            bool signal = false;
            double timeRef = 0;
            
            int nEntries = tSignal->GetEntries();
            for (int l = 0; l < nEntries; l++) {
                    tSignal->GetEntry(l);
                    if (ft < timeRef) {
                        std::cout << "erreur : le temps n'est pas lu dans l'ordre croissant" << std::endl;
                        return 0;
                    }
                double tStep = ft-timeRef;
                //std::cout << "tStep = " << tStep << std::endl;
                timeRef = ft;
                if (abs(fct)>threshold) {
                    if (!signal) {nEvents++; signal = true;}
                    nIons+=abs(fci)/tStep;
                    nElectrons+=abs(fce)/tStep;
                    nTotal+=abs(fct)/tStep;
                }
                else {
                    if (signal) {
                        tCharges->Fill();
                        nIons = 0; nElectrons = 0; nTotal = 0;
                        signal = false;
                    }
                    
                }
                /*
                if (abs(chargeNumberIon) > 0.0005) hIon->Fill(abs(chargeNumberIon));
                if (abs(chargeNumberElectron) > 10) hElectron->Fill(abs(chargeNumberElectron));
                if (abs(chargeNumberTotal) > 10) hTotal->Fill(abs(chargeNumberTotal));
                
                std::cout << "chargeNumberIon " << chargeNumberIon << std::endl;
                std::cout << "chargeNumberElectron " << chargeNumberElectron << std::endl;
                std::cout << "chargeNumberTotal " << chargeNumberTotal << std::endl;
                 */
                
            }
            /*
            c2->cd(1);
            hIon->Draw();
            hIon->Write();
            c2->cd(2);
            hElectron->Draw();
            hElectron->Write();
            c2->cd(3);
            hTotal->Draw();
            hTotal->Write();
             */
            tCharges->Write();
            
            //c2->SaveAs(Form("Figures/CurrentsV%d.pdf", j+2));
            
            std::cout << nEvents << " events founds " << std::endl;
            

        }
        fAmp->Close();
     //}
    
     
    //cvf->SaveAs("Figures/ElectronCurrents.pdf");
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    return 0;
}


