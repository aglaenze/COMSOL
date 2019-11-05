#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "_Utils.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/ViewFEMesh.hh"
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
    
    if (argc < 2) {
        std::cout << "Please enter HVmesh like this: ./avalanche $hvMesh " << std::endl;
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
    //const double radius = 0.0004;   // cm
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
    
    // To look at the avalanche
    ViewDrift* driftView = new ViewDrift();
    ViewDrift* driftView2 = new ViewDrift();
    driftView->SetArea(depth/2. - 10*pitch, depth/2. - 10*pitch, 0, depth/2. + 10*pitch, depth/2. + 10*pitch, damp*3);
    driftView2->SetArea(depth/2. - 5*pitch, depth/2. - 5*pitch, 0, depth/2. + 5*pitch, depth/2. + 5*pitch, damp*2);
    
    // Create an avalanche object
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(&sensor);
    AvalancheMC* drift = new AvalancheMC();
    drift->SetSensor(&sensor);
    //drift->SetDistanceSteps(2.e-4);
    aval->EnablePlotting(driftView);
    aval->EnablePlotting(driftView2);
    drift->EnablePlotting(driftView);
    drift->EnablePlotting(driftView2);
    
    const int nEvents = 1;
    
    for (unsigned int i = 0; i < nEvents; ++i) {
        // Initial coordinates of the electron.
        double x0 = width/2. + RndmUniform() * pitch;
        //double y0 = RndmUniform() * depth;
        double y0 = depth/2. + RndmUniform() * pitch;
        //double z0 = 2*damp + (ddrift-2*damp)*RndmUniform();
        double z0 = 2*damp;
        double t0 = 0;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        int ne2 = 0, ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        if (ne2==1) i--;
        std::cout << "\nAvalanche size = " << ne2 << std::endl;
        if (ne2==1) {i--; continue;}
        //std::cout << "where is the problem??" << std::endl;
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        double xi1, yi1, zi1, ti1;
        double xi2, yi2, zi2, ti2;
        int status;
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
                                      xe2, ye2, ze2, te2, e2, status);
            drift->DriftIon(xe1, ye1, ze1, te1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
        }
    }
    
    const bool plotDrift = true;
    if (plotDrift) {
        TCanvas* c1 = new TCanvas();
        driftView->SetCanvas(c1);
        driftView->Plot();
        c1->SaveAs(Form("Figures/avalanche_%s.pdf", gasName.c_str()));
    }
  
        
    // Set up the object for FE mesh visualization.
    ViewFEMesh* vFE = new ViewFEMesh();
    vFE->SetArea(depth/2. - 5*pitch, depth/2. - 5*pitch, 0, depth/2. + 5*pitch, depth/2. + 5*pitch, damp*2);
    vFE->SetComponent(fm);
    vFE->SetViewDrift(driftView2);
    vFE->SetPlane(0, -1, 0, width/2., depth/2., damp);
    vFE->SetFillMesh(true);
    const bool plotDrift2 = true;
    if (plotDrift2) {
        TCanvas* c2 = new TCanvas();
        vFE->SetCanvas(c2);
        //for (unsigned int i = 0; i < nMaterials; ++i) vFE->SetColor(i, kGray+i);
        //vFE->SetColor(0, kWhite);
        //vFE->SetColor(0, kYellow + 3);
        vFE->SetColor(1, kGray);
        //vFE->SetColor(2, kYellow + 3);
        vFE->EnableAxes();
        vFE->SetXaxisTitle("x (cm)");
        vFE->SetYaxisTitle("z (cm)");
        std::cout << "Plotting..." << std::endl;
        vFE->Plot();
        c2->SaveAs(Form("Figures/avalanche2d_%s.pdf", gasName.c_str()));
    }
    
    time_t t1 = time(NULL);
    
    //std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    app.Run(true);
    
    //return 0;
}


