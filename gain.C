#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "_Utils.C"
#include "parameters.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
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
    const int modelNum = 10;
    //____________________
    
    time_t t0 = time(NULL);
    
    TString fOutputName;
    
    // Make a gas medium.
    MediumMagboltz* gas = InitiateGas(gasName);
    // Load field map
    ComponentComsol* fm;
    
    int hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvDrift = 0, saveNum = 0;
    if (modelNum == 1) {
        if (argc < 4) {
            std::cout << "Please enter HVmesh like this: ./gain $hvMesh $saveNum " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDrift = atoi(argv[2]);
        saveNum = atoi(argv[3]);
        fm = InitiateField(modelNum, hvMesh, hvDrift, gas);
        fOutputName = Form("Gain-%d-%d-model%d-%d.pdf", hvMesh, hvDrift, modelNum, saveNum);
    }
    else if (modelNum > 6 && modelNum < 10) {
        if (argc < 5) {
            std::cout << "Please enter HVmesh like this: ./gain $hvDmDown $hvDmUp $saveNum " << std::endl;
            return 0;
        }
        hvDmDown = atoi(argv[1]);
        hvDmUp = atoi(argv[2]);
        hvDrift = atoi(argv[3]);
        saveNum = atoi(argv[4]);
        fm = InitiateField(modelNum, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("Gain-%d-%d-%d-model%d-%d.pdf", hvDmDown, hvDmUp, hvDrift, modelNum, saveNum);
    }
    else if (modelNum == 10) {
        if (argc < 6) {
            std::cout << "Please enter HVmesh like this: ./gain $hvMesh $hvDmDown $hvDmUp $saveNum " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDmDown = atoi(argv[2]);
        hvDmUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        saveNum = atoi(argv[5]);
        fm = InitiateField(modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("Gain-%d-%d-%d-%d-model%d-%d.pdf", hvMesh, hvDmDown, hvDmUp, hvDrift, modelNum, saveNum);
    }
    else {std::cout << "Wrong model number" << std::endl; return 0;}

    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
        //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
    //std::cout << damp << " " << width << " " << depth << " " << ddrift << std::endl;
    //return 0;

    
    // Make a sensor.
    Sensor sensor;
    sensor.AddComponent(fm);
    sensor.SetArea(0, 0, 0, width, depth, ddrift);
    //sensor.SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);

    // Create an avalanche object
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(&sensor);

    // Create ROOT histograms of the signal and a file in which to store them.
    //const int nBins = 50000;    //12000
    int nBins;
    //if (modelNum == 1) nBins = int(0.2*TMath::Exp(0.0352*hvMesh));
    if (modelNum == 1) nBins = 50000;
    //else if (modelNum > 6) nBins = int(0.2*TMath::Exp(0.0352*hvDmUp));
    else if (modelNum > 6) nBins = 50000;
    TH1F* hElectrons = new TH1F("hElectrons", "Number of secondary electrons", nBins, 0, 4*nBins);
    //TH1::StatOverflows(true);
    hElectrons->SetXTitle("# secondary electrons");
    hElectrons->SetYTitle("# counts");
    
    // Write the histograms to the TFile.
    char* name;
    if (modelNum == 1) name = Form("rootFiles/%s/model%d/gain-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDrift, saveNum);
    else if (modelNum > 6 && modelNum < 10) name = Form("rootFiles/%s/model%d/gain-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift, saveNum);
    else if (modelNum == 10 ) name = Form("rootFiles/%s/model%d/gain-%d-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, saveNum);
    TFile* f = new TFile(name, "RECREATE");
    
    const int nEvents = 10000;
    //const int nEvents = 5000;
    int division = int(nEvents/20);
    
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (i % 100 == 0) {
            std::cout << "\n\n\n\n" << i << "/" << nEvents << std::endl;
            time_t t = time(NULL);
            PrintTime(t0, t);            
        }
        if (i % division == 0) hElectrons->Write("", TObject::kOverwrite);
        // Initial coordinates of the photon.
        double x0 = width/2. + RndmUniform() * pitch;
        //double y0 = RndmUniform() * depth;
        double y0 = depth/2. + RndmUniform() * pitch;
        double z0 = damp + radius + (ddrift-damp-radius)*RndmUniform();
        //double z0 = damp + 3*radius;
        double t0 = 0;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        int ne2 = 0, ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        std::cout << "Avalanche size = " << ne2 << std::endl;
        if (ne2 < 2) continue;
        if (modelNum == 1) {
            hElectrons->Fill(ne2);    // ok for modelNum < 4
            continue;
        }
        // need to look at the end points of the electron in the avalanche for modelNum >= 4
        // if electrons are stopped above the last MM, they won't count in the gain
        // need to select the winners
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        int status;
        int nWinners = 0;
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
            //std::cout << "departure of the electron in x y z : " << xe1 << " " << ye1 << " " <<  ze1 << std::endl;
            //std::cout << "arrival of the electron in x y z : " << xe2 << " " << ye2 << " " <<  ze2 << std::endl;
            if (ze2 < 0.01) nWinners++;
        }
        std::cout << "nWinners = " << nWinners << " / " << ne2 << std::endl;
        if (nWinners > 0) hElectrons->Fill(nWinners);
    }
    f->Close();
    
    const bool drawSpectrum = false;
    if (drawSpectrum) {
        TCanvas* c = new TCanvas();
        c->cd();
        hElectrons->SetFillColor(kBlue + 2);
        hElectrons->SetLineColor(kBlue + 2);
        hElectrons->Draw();
        c->SaveAs(fOutputName);
    }
    
    
    time_t t1 = time(NULL);
    
    std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
}



