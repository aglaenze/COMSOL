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
    std::string gasName = "Ar"; // Ar or Ne
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
    
    // Set up detector geometry
    //GeometrySimple* geo = new GeometrySimple();
    const double pitch = 0.0025;    // cm
    const double damp = 0.0128;
    const double ddrift = 0.5;      // cm
    const double radius = 0.0004;   // cm
    const int periodicityNum = 5000;
    double width = periodicityNum * pitch;
    double depth = periodicityNum * pitch;
    
    
    // Make a gas medium.
    MediumMagboltz* gas = new MediumMagboltz();
    double rPenning;
    const double lambdaPenning = 0.;    // parameter for sampling the distance of the Penning electron with respect to the excitation
    const std::string path = getenv("GARFIELD_HOME");
    if (gasName=="Ar") {
        gas->SetComposition("Ar", 95., "iC4H10", 5.);
        rPenning = 0.45;
        gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
        gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
    }
    else if (gasName=="Ne") {
        gas->SetComposition("Ne", 90., "CF4", 10.);
        rPenning = 0.6;
        gas->EnablePenningTransfer(rPenning, lambdaPenning, "ne");
        gas->LoadIonMobility(path + "/Data/IonMobility_Ne+_Ne.txt");
    }
    else {std::cout << "What gas??" << std::endl; return 0;}
    //return 0;
    gas->SetTemperature(293.15);
    gas->SetPressure(AtmosphericPressure);
    gas->EnableDrift();

    //gas->EnablePenningTransfer(rPenning, lambdaPenning);
    
    // Load the field map.
    std::string dataFolder = Form("COMSOL_data/model%d/", modelNum);
    std::string dataFile = dataFolder + Form("ewfield_%dV.txt", hvMesh);
    // Load the field map.
    ComponentComsol* fm = new ComponentComsol();
    fm->Initialise(dataFolder+"mesh.mphtxt", dataFolder+"dielectrics.dat", dataFile);
    fm->PrintMaterials();
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
    fm->PrintRange();

    
    // Associate the gas with the corresponding field map material.
    /*
     const unsigned int nMaterials = fm->GetNumberOfMaterials();
     for (unsigned int i = 0; i < nMaterials; ++i) {
     const double eps = fm->GetPermittivity(i);
     if (eps == 1.) fm->SetMedium(i, gas);
     }
     */
    fm->SetMedium(0, gas);
    
    fm->PrintMaterials();
    
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
    int nBins = int(0.4*TMath::Exp(0.0352*hvMesh));
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



