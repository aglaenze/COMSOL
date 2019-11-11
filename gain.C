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
    std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    //____________________
    
    time_t t0 = time(NULL);
    
    /*
    if (argc < 3) {
        std::cout << "Please enter HVmesh like this: ./gain $hvMesh $gain$i " << std::endl;
        return 0;
    }
     */
    if (argc < 2) {
        std::cout << "Please enter HVmesh like this: ./gain $hvMesh " << std::endl;
        return 0;
    }
    
    const int hvMesh = atoi(argv[1]);
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();

    // Make a gas medium.
    MediumMagboltz* gas = InitiateGas(gasName);
    // Load field map
    ComponentComsol* fm = InitiateField(modelNum, hvMesh, gas);
    
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
    int nBins = int(0.05*TMath::Exp(0.0352*hvMesh));
    TH1F* hElectrons = new TH1F("hElectrons", "Number of secondary electrons", int(nBins/4.), 0, nBins);
    //TH1::StatOverflows(true);
    hElectrons->SetXTitle("# secondary electrons");
    hElectrons->SetYTitle("# counts");
    
    // Write the histograms to the TFile.
    //const char* name = Form("rootFiles/%s/model%d/gain_%dV_%s.root", gasName.c_str(), modelNum, hvMesh, argv[2]);
    const char* name = Form("rootFiles/%s/model%d/gain_%dV.root", gasName.c_str(), modelNum, hvMesh);
    TFile* f = new TFile(name, "RECREATE");
    
    const int nEvents = 12000;
    int division = int(nEvents/20);
    
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (i % 100 == 0) std::cout << "\n\n\n\n" << i << "/" << nEvents << "\n";
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
        hElectrons->Fill(ne2);
        if (i % division == 0) hElectrons->Write("", TObject::kOverwrite);
    }
    f->Close();
    
    const bool drawSpectrum = false;
    if (drawSpectrum) {
        TCanvas* c = new TCanvas();
        c->cd();
        hElectrons->SetFillColor(kBlue + 2);
        hElectrons->SetLineColor(kBlue + 2);
        hElectrons->Draw();
        c->SaveAs(Form("Gain_%dV_model%d.pdf", hvMesh, modelNum));
    }
    
    
    time_t t1 = time(NULL);
    
    std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
}



