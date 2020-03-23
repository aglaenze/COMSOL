/*
Note that the spectrum shape only depends on the detector geometry
 Photons create electrons with energy = 0
*/

#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "_Utils.C"
#include "parameters.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {
    
    //______________________
    // variables
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
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
    int saveNum;
    if (modelNum == 1) {
        if (argc != 4) {
            std::cout << "Please enter HVmesh like this: ./feSignal $hvMesh $hvDrift $saveNum" << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDrift = atoi(argv[2]);
        saveNum = atoi(argv[3]);
        fm = InitiateField(modelNum, hvMesh, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/feSignal-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDrift, saveNum);
    }
    else if (modelNum >= 2 && modelNum < 5) {
        if (argc != 5) {
            std::cout << "Please enter HVmesh like this: ./feSignal $hvDmDown $hvDmUp $hvDrift $saveNum" << std::endl;
            return 0;
        }
        hvDmDown = atoi(argv[1]);
        hvDmUp = atoi(argv[2]);
        hvDrift = atoi(argv[3]);
        saveNum = atoi(argv[4]);
        fm = InitiateField(modelNum, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/feSignal-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift, saveNum);
    }
    else if (modelNum >= 5 && modelNum < 8) {
        if (argc != 6) {
            std::cout << "Please enter HVmesh like this: ./feSignal $hvMesh $hvDmDown $hvDmUp $hvDrift $saveNum" << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDmDown = atoi(argv[2]);
        hvDmUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        saveNum = atoi(argv[5]);
        fm = InitiateField(modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/feSignal-%d-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, saveNum);
    }
    else if (modelNum >= 8 && modelNum < 10) {
        if (argc != 6) {
            std::cout << "Please enter HVmesh like this: ./feSignal $hvMesh $hvGemDown $hvGemUp $hvDrift $saveNum" << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvGemDown = atoi(argv[2]);
        hvGemUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        saveNum = atoi(argv[5]);
        fm = InitiateField(modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/feSignal-%d-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift, saveNum);
    }
    else {std::cout << "Wrong model number" << std::endl; return 0;}
    
    // for testing without overwriting good files
    //fOutputName = Form("rootFiles/%s/model%d/test-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift);
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    int electrodeNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth, electrodeNum);


    // Set up detector geometry
    GeometrySimple* geo = new GeometrySimple();
    // Set up the bulk
    SolidBox bulk = SolidBox(width/2., depth/2., ddrift/2., width/2., depth/2., ddrift/2.);
    geo->AddSolid(&bulk, gas);
    fm->SetGeometry(geo);
    geo->PrintSolids();

    // Make a sensor
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
    sensor->SetArea(0, 0, 0, width, depth, ddrift);
    for (int k = 0; k < electrodeNum; k++) sensor->AddElectrode(fm, Form("V%d", k+2));
    
    // Use Heed for simulating the photon absorption.
    TrackHeed track;
    track.SetSensor(sensor);
    track.EnableElectricField();
    
    const int nEvents = 100;
    
    // Create ROOT histograms of the signal and a file in which to store them.
    TFile* f = new TFile(fOutputName, "RECREATE");
    //TFile* f = new TFile("rootFiles/test.root", "RECREATE");
    
    TTree *tGain = new TTree("tGain","Gain");
    Int_t nWinners = 0;
    tGain->Branch("electronNumber", &nWinners, "electronNumber/I");

    // Set the signal binning.
    //const int nTimeBins = 10000;
    const double tStep = 0.1;   //ns
    const double rate = 6.e7;               // number of events per s (note that t units are ns here, we'll need a conversion factor)
    const double timespace = 1./rate*1.e9;    // in ns
    const double ionDelay = 1.e3;    // time to collect all ions at the drift electrode ~1ms
    const double tStart =  0.;
    //const double tEnd = int(nEvents * timespace + ionDelay);
    const double tEnd = int(nEvents * timespace + ionDelay);
    //const double tStep = (tEnd - tStart) / nTimeBins;
    const int nTimeBins = (tEnd - tStart)/tStep;
    sensor->SetTimeWindow(tStart, tStep, nTimeBins);
    
    // To look at the avalanche
    ViewDrift* driftView = new ViewDrift();
    //driftView->SetArea(390*pitch, 390*pitch, 0, 410*pitch, 410*pitch, damp+6*pitch);
    driftView->SetArea(0, 0, 0, width, depth, ddrift);
    
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
            std::cout << "\n\n\n\n" << i << "/" << nEvents << "\n\n" << std::endl;
            time_t t = time(NULL);
            PrintTime(t0, t);
            tGain->Write("", TObject::kOverwrite);
        }
        // Initial coordinates of the photon.
        double x0 = width/2.;
        double y0 = depth/2.;
        double z0 = ddrift;
        double t0 = (i+RndmUniform() ) * timespace;
        // Sample the photon energy, using the relative intensities according to XDB.
        const double r = 167. * RndmUniform();
        const double egamma = r < 100. ? 5898.8 : r < 150. ? 5887.6 : 6490.4;
        /* // translation
         double egamma;
         if ( r<100) egamma = 5898.8;
         else if ( r<150 ) egamma = 5887.6;
         else egamma = 6490.4;
         */
        int ne = 0;
        double dx0 = -1+RndmUniform()*2;
        double dy0 = -1+RndmUniform()*2;
        track.TransportPhoton(x0, y0, z0, t0, egamma, dx0, dy0, -1, ne);
        //track.TransportPhoton(x0, y0, z0, t0, egamma, 0, 0, -1, ne);
        if (ne < 2) {i--; continue;}
        std::cout << "ne = " << ne << std::endl;    // number of primaries
        //continue;
        nWinners = 0;
        double x, y, z, t, e, dx, dy, dz;
        for (int j = 0; j < ne; j++) {
            track.GetElectron(j, x, y, z, t, e, dx, dy, dz);
            //std::cout << x << " " << y << " " << z << " " << t << " " << e << " " << std::endl;
            aval->AvalancheElectron(x, y, z, t, e, dx, dy, dz);
            int ne2 = 0, ni = 0;
            aval->GetAvalancheSize(ne2, ni);
            // std::cout << "\nAvalanche size = " << ne2 << std::endl;
            //if (ne2 < 4) {j--; continue;}
            /*
            if (modelNum == 1) { nWinners = ne2;
                //hElectrons->Fill(ne2);    // ok for modelNum < 4
                //continue;
            }
             */
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
                if (ze2 < 0.01) nWinners++;
            }
        }
        tGain->Fill();
        
    }
    tGain->Write();

    for (int k = 0; k < electrodeNum; k++) {
        TTree *tSignal = new TTree(Form("tSignal_V%d",k+2),"Currents");
        Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
        tSignal->Branch("time", &ft, "time/D");
        tSignal->Branch("totalCurrent", &fct, "totalCurrent/D");
        tSignal->Branch("electronCurrent", &fce, "electronCurrent/D");
        tSignal->Branch("ionCurrent", &fci, "ionCurrent/D");

        
        for (int j = 0; j < nTimeBins; j++) {
            ft = j * tStep;
            //std::cout << ft << std::endl;
            fct = sensor->GetSignal(Form("V%d", k+2), j) / ElementaryCharge;
            fce = sensor->GetElectronSignal(Form("V%d", k+2), j) / ElementaryCharge;
            fci = sensor->GetIonSignal(Form("V%d", k+2), j) / ElementaryCharge;
            tSignal->Fill();
        }
        tSignal->Write();
    }

    f->Close();
    
    // Plot the signal.
    const bool plotSignal = false;
    if (plotSignal) {
        TCanvas* cSignal = new TCanvas("feSignal", "Fe Signal");
        ViewSignal* vSignal = new ViewSignal();
        gPad->SetLeftMargin(0.15);
        vSignal->SetSensor(sensor);
        vSignal->SetCanvas(cSignal);
        //void PlotSignal(const std::string& label, const bool total = true, const bool electron = false, const bool ion = false);
        vSignal->SetRangeX(0, tEnd);
        vSignal->PlotSignal("V2", true, true, true);    // signal de la mesh
        cSignal->SaveAs("Figures/FeSignal.pdf");
    }
    
    time_t t1 = time(NULL);
    
    std::cout << "\n" << nEvents << " events of photons simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
}


