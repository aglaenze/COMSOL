#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "_Utils.C"
#include "parameters.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 8;
    //____________________
    
    time_t t0 = time(NULL);
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    TString fOutputName;
    
    // Make a gas medium.
    MediumMagboltz* gas = InitiateGas(gasName);
    // Load field map
    ComponentComsol* fm;
    
    int hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvGemDown = 0, hvGemUp = 0, hvDrift = 0;
    int electrodeNum = 0;
    if (modelNum == 1) {
        if (argc != 3) {
            std::cout << "Please enter HVmesh like this: ./signal $hvMesh $hvDrift " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDrift = atoi(argv[2]);
        fm = InitiateField(modelNum, hvMesh, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDrift);
        electrodeNum = 3;
    }
    else if (modelNum >= 2 && modelNum < 5) {
        if (argc != 4) {
            std::cout << "Please enter HVmesh like this: ./signal $hvDmDown $hvDmUp $hvDrift " << std::endl;
            return 0;
        }
        hvDmDown = atoi(argv[1]);
        hvDmUp = atoi(argv[2]);
        hvDrift = atoi(argv[3]);
        fm = InitiateField(modelNum, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d.root", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift);
        electrodeNum = 4;
    }
    else if (modelNum >= 5 && modelNum < 8) {
        if (argc != 5) {
            std::cout << "Please enter HVmesh like this: ./signal $hvMesh $hvDmDown $hvDmUp $hvDrift " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDmDown = atoi(argv[2]);
        hvDmUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        fm = InitiateField(modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift);
        electrodeNum = 5;
    }
    else if (modelNum >= 8 && modelNum < 10) {
        if (argc != 5) {
            std::cout << "Please enter HVmesh like this: ./signal $hvMesh $hvGemDown $hvGemUp $hvDrift " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvGemDown = atoi(argv[2]);
        hvGemUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        fm = InitiateField(modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift);
        electrodeNum = 5;
    }
    else {std::cout << "Wrong model number" << std::endl; return 0;}
    
        //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
    //std::cout << damp << " " << width << " " << depth << " " << ddrift << std::endl;
    //return 0;

    
    // Make a sensor.
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
    sensor->SetArea(0, 0, 0, width, depth, ddrift);
    //sensor->SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
    for (int k = 0; k < electrodeNum; k++) sensor->AddElectrode(fm, Form("V%d", k+1));
    
    // Set the signal binning.
    //const int nTimeBins = 10000;
    const double tStep = 0.1;   //ns
    const int nEvents = 700;
    const double rate = 1.e7;               // number of events per s (note that t units are ns here)
    const double timespace = 1./rate*1.e9;    // in ns
    const double ionDelay = 1.e3;    // time to collect all ions at the drift electrode ~1ms
    const double tStart =  0.;
    const double tEnd = int(nEvents * timespace + ionDelay);
    //const double tStep = (tEnd - tStart) / nTimeBins;
    const int nTimeBins = (tEnd - tStart)/tStep;
    sensor->SetTimeWindow(tStart, tStep, nTimeBins);
    
    // To look at the avalanche
    ViewDrift* driftView = new ViewDrift();
    driftView->SetArea(390*pitch, 390*pitch, 0, 410*pitch, 410*pitch, damp+6*pitch);
    
    // Create an avalanche object
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    aval->EnablePlotting(driftView);
    aval->EnableSignalCalculation();
    
    AvalancheMC* drift = new AvalancheMC();
    drift->SetSensor(sensor);
    drift->SetDistanceSteps(2.e-4);
    drift->EnablePlotting(driftView);
    drift->EnableSignalCalculation();
    
    int division = int(nEvents/10);
    
    
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (i % division == 0) {
            std::cout << "\n\n\n\n" << i << "/" << nEvents << std::endl;
            time_t t = time(NULL);
            PrintTime(t0, t);            
        }
        // Initial coordinates of the photon.
        double x0 = width/2. + RndmUniform() * pitch;
        //double y0 = RndmUniform() * depth;
        double y0 = depth/2. + RndmUniform() * pitch;
        double z0 = damp + 2*radius + (ddrift-damp-2*radius)*RndmUniform();
        double t0 = ( i + RndmUniform() )* timespace;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        int ne2 = 0, ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        std::cout << "\nAvalanche size = " << ne2 << std::endl;
        if (ne2 < 4) {i--; continue;}
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        double xi1, yi1, zi1, ti1;
        double xi2, yi2, zi2, ti2;
        int status;
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
            drift->DriftIon(xe1, ye1, ze1, te1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
        }
    }
    
    const bool extract = true;
    if (extract) {

        // Create ROOT histograms of the signal and a file in which to store them.
        TFile* f = new TFile(fOutputName, "RECREATE");


        for (int k = 0; k < electrodeNum; k++) {
            TTree *tSignal = new TTree(Form("tSignal_V%d",k+1),"Currents");
            Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
            tSignal->Branch("time", &ft, "time/D");
            tSignal->Branch("totalCurrent", &fct, "totalCurrent/D");
            tSignal->Branch("electronCurrent", &fce, "electronCurrent/D");
            tSignal->Branch("ionCurrent", &fci, "ionCurrent/D");
            
            /*
            //TList* fOutputList = new TList();          // this is a list which will contain all  histograms
            //fOutputList->SetOwner(kTRUE);
            TH1F* hS = new TH1F("htotalCurrent", "Total signal (mesh)", nTimeBins, 0, tEnd);
            TH1F* hSE = new TH1F("helectronCurrent", "Electron signal (mesh)", nTimeBins, 0, tEnd);
            TH1F* hSI = new TH1F("hionCurrent", "Ion signal (mesh)", nTimeBins, 0, tEnd);
            
            tSignal->Branch("htotalCurrent", &hS, "hTotalCurrent");
            tSignal->Branch("helectronCurrent", &hSE, "hElectronCurrent");
            tSignal->Branch("hionCurrent", &hSI, "hIonCurrent");
            
            //fOutputList->Add(hS);
            //fOutputList->Add(hSE);
            //fOutputList->Add(hSI);
             */
             
            for (int j = 0; j < nTimeBins; j++) {
                ft = j * tStep;
                //std::cout << ft << std::endl;
                fct = sensor->GetSignal(Form("V%d", k+1), j) / ElementaryCharge;
                fce = sensor->GetElectronSignal(Form("V%d", k+1), j) / ElementaryCharge;
                fci = sensor->GetIonSignal(Form("V%d", k+1), j) / ElementaryCharge;
                /*
                hS->Fill(j * tStep, fct);
                hSE->Fill(j * tStep, fce);
                hSI->Fill(j * tStep, fci);
                 */
                tSignal->Fill();
            }
            tSignal->Write();
            //fOutputList->Write();
        }
        
        
        /*
         // Extract the calculated signal.
         double bscale = tEnd / nTimeBins;  // time per bin
         double sum = 0.;  // to keep a running sum of the integrated signal
         double sumE = 0.;
         double sumI = 0.;
         double sumd = 0.;  // drift electrode
         double sumEd = 0.;
         double sumId = 0.;
         */
          
        /*
         TH1F* hS = new TH1F("hhMesh", "Total signal (mesh)", nTimeBins, 0, tEnd);
         TH1F* hInt = new TH1F("hIntMesh", "Integrated signal (mesh)", nTimeBins, 0, tEnd);
        for (int i = 0; i < nTimeBins; i++) {
            //double wt = sensor->GetSignal("V2", i) / ElementaryCharge;
            //sum += wt;
            //hS->Fill(i * bscale, wt);
            hS->Fill(1);
            hInt->Fill(1);
        }
        
        TList* fOutputList = new TList();          // this is a list which will contain all  histograms
        fOutputList->Add(hS);
        fOutputList->Add(hInt);
        fOutputList->Write();
        //return 0;
         */
        
        
        /*
        TH1F* hS = new TH1F("hhMesh", "Total signal (mesh)", nTimeBins, 0, tEnd);
        TH1F* hInt = new TH1F("hIntMesh", "Integrated signal (mesh)", nTimeBins, 0, tEnd);
        TH1F* hSE = new TH1F("hhEMesh", "Electron signal (mesh)", nTimeBins, 0, tEnd);
        TH1F* hIntE = new TH1F("hIntEMesh", "Integrated electron signal (mesh)", nTimeBins, 0, tEnd);
        TH1F* hSI = new TH1F("hhIMesh", "Ion signal (mesh)", nTimeBins, 0, tEnd);
        TH1F* hIntI = new TH1F("hIntIMesh", "Integrated ion signal (mesh)", nTimeBins, 0, tEnd);

        TH1F* hSd = new TH1F("hhDrift", "Total signal (drift)", nTimeBins, 0, tEnd);
        TH1F* hIntd = new TH1F("hIntDrift", "Integrated signal (drift)", nTimeBins, 0, tEnd);
        TH1F* hSEd = new TH1F("hhEDrift", "Electron signal (drift)", nTimeBins, 0, tEnd);
        TH1F* hIntEd = new TH1F("hIntEDrift", "Integrated electron signal (drift)", nTimeBins, 0, tEnd);
        TH1F* hSId = new TH1F("hhIDrift", "Ion signal (drift)", nTimeBins, 0, tEnd);
        TH1F* hIntId = new TH1F("hIntIDrift", "Integrated ion signal (drift)", nTimeBins, 0, tEnd);

        // Fill the histograms with the signals.
        //  Note that the signals will be in C/(ns*binWidth), and we will divide by e
        // to give a signal in e/(ns*binWidth).
        //  The total signal is then the integral over all bins multiplied by the bin
        // width in ns.
        for (int i = 0; i < nTimeBins; i++) {
            double wt = sensor->GetSignal("V2", i) / ElementaryCharge;
            sum += wt;
            hS->Fill(i * tStep, wt);
            hInt->Fill(i * tStep, sum);
            double wtE = sensor->GetElectronSignal("V2", i) / ElementaryCharge;
            sumE += wtE;
            hSE->Fill(i * tStep, wtE);
            hIntE->Fill(i * tStep, sumE);
            double wtI = sensor->GetIonSignal("V2", i) / ElementaryCharge;
            sumI += wtI;
            hSI->Fill(i * tStep, wtI);
            hIntI->Fill(i * tStep, sumI);
            
            double wtd = sensor->GetSignal("V3", i) / ElementaryCharge;
            sumd += wtd;
            hSd->Fill(i * tStep, wtd);
            hIntd->Fill(i * tStep, sumd);
            double wtEd = sensor->GetElectronSignal("V3", i) / ElementaryCharge;
            sumEd += wtEd;
            hSEd->Fill(i * tStep, wtEd);
            hIntEd->Fill(i * tStep, sumEd);
            double wtId = sensor->GetIonSignal("V3", i) / ElementaryCharge;
            sumId += wtId;
            hSId->Fill(i * tStep, wtId);
            hIntId->Fill(i * tStep, sumId);
        }
        //std::cout << "IBF = " << sumd*100./sum - 100./gain << "%" << std::endl;
        std::cout << "IBF = " << sumd*100./sum << "%" << std::endl;
        
        hS->SetOption("hist");
        hInt->SetOption("hist");
        hSE->SetOption("hist");
        hIntE->SetOption("hist");
        hSI->SetOption("hist");
        hIntI->SetOption("hist");
        hSd->SetOption("hist");   // drift electrode
        hIntd->SetOption("hist");
        hSEd->SetOption("hist");
        hIntEd->SetOption("hist");
        hSId->SetOption("hist");
        hIntId->SetOption("hist");
        
        hS->SetXTitle("time (ns)");
        hS->SetYTitle("# charges / unit time");
        hSE->SetXTitle("time (ns)");
        hSE->SetYTitle("# charges / unit time");
        hSI->SetXTitle("time (ns)");
        hSI->SetYTitle("# charges / unit time");
        hSd->SetXTitle("time (ns)");
        hSd->SetYTitle("# charges / unit time");
        hSEd->SetXTitle("time (ns)");
        hSEd->SetYTitle("# charges / unit time");
        hSId->SetXTitle("time (ns)");
        hSId->SetYTitle("# charges / unit time");
        
        
        // Write the histograms to the TFile.
        hS->Write();    // mesh electrode
        hInt->Write();
        hSE->Write();
        hIntE->Write();
        hSI->Write();
        hIntI->Write();

        hSd->Write();   // drift electrode
        hIntd->Write();
        hSEd->Write();
        hIntEd->Write();
        hSId->Write();
        hIntId->Write();
         */
        

        //f->Write();
        f->Close();
    }
    
    // Plot the signal.
    const bool plotSignal = false;
    if (plotSignal) {
        TCanvas* cSignal = new TCanvas("signal", "Signal");
        ViewSignal* vSignal = new ViewSignal();
        gPad->SetLeftMargin(0.15);
        vSignal->SetSensor(sensor);
        vSignal->SetCanvas(cSignal);
        //void PlotSignal(const std::string& label, const bool total = true, const bool electron = false, const bool ion = false);
        vSignal->SetRangeX(0, nEvents*timespace+200);
        vSignal->PlotSignal("V2", true, true, true);    // signal de la mesh
        cSignal->SaveAs("Figures/Signal.pdf");
    }
    
    
    time_t t1 = time(NULL);
    
    std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
}


