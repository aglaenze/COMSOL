#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "_Utils.C"
#include "parameters.C"

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
    std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 6;
    //____________________
    
    time_t t0 = time(NULL);
    
    if (argc < 2) {
        std::cout << "Please enter HVmesh like this: ./avalanche $hvMesh " << std::endl;
        return 0;
    }
    const int hvMesh = atoi(argv[1]);
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    //Load geometry parameters
        //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
    
    // Make a gas medium.
    MediumMagboltz* gas = InitiateGas(gasName);
    // Load field map
    ComponentComsol* fm = InitiateField(modelNum, hvMesh, gas);
    
    // Make a sensor.
    Sensor sensor;
    sensor.AddComponent(fm);
    sensor.SetArea(0, 0, 0, width, depth, ddrift);
    //sensor.SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
    
    // To look at the avalanche
    ViewDrift* driftView = new ViewDrift();
    ViewDrift* driftView2 = new ViewDrift();
    driftView->SetArea(-10*pitch, -10*pitch, 0, 10*pitch, 10*pitch, damp*3);
    driftView2->SetArea(-5*pitch, -5*pitch, 0, 5*pitch, 5*pitch, damp*2);
    
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
    
    int j = 0;
    for (unsigned int i = 0; i < nEvents; ++i) {
        j++;
        if (j==10) break;
        //std::cout << "hello " << i << "\n\n" << std::endl;
        // Initial coordinates of the electron.
        double x0 = RndmUniform() * pitch;
        double y0 = RndmUniform() * depth;
        //double z0 = 2*damp + (ddrift-2*damp)*RndmUniform();
        double z0 = 0.3;
        double t0 = 0;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        int ne2 = 0, ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        std::cout << "\nAvalanche size = " << ne2 << std::endl;
        if (ne2 < 100) {i--; continue;}
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        double xi1, yi1, zi1, ti1;
        double xi2, yi2, zi2, ti2;
        int status;
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
            /*
            if (ze2 < 0.012){
            std::cout << "departure of the electron in x y z : " << xe1 << " " << ye1 << " " <<  ze1 << std::endl;
            std::cout << "arrival of the electron in x y z : " << xe2 << " " << ye2 << " " <<  ze2 << std::endl;
            }
             */
            drift->DriftIon(xe1, ye1, ze1, te1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
        }
    }
    
    return 0;
    
    const bool plotDrift = false;
    if (plotDrift) {
        TCanvas* c1 = new TCanvas();
        driftView->SetCanvas(c1);
        driftView->Plot();
        c1->SaveAs(Form("Figures/avalanche_%s.pdf", gasName.c_str()));
    }
  
        
    // Set up the object for FE mesh visualization.
    ViewFEMesh* vFE = new ViewFEMesh();
    //vFE->SetArea(-5*pitch, -5*pitch, 0,  5*pitch, 5*pitch, damp*2);
    //vFE->SetArea(-5*pitch, -5*pitch, 0.2-2*pitch,  5*pitch, 5*pitch, 0.2+8*pitch);
    vFE->SetArea(-0.15, -0.15, 0,  0.15, 0.15, 0.3);
    vFE->SetComponent(fm);
    vFE->SetViewDrift(driftView2);
    //vFE->SetPlane(0, -1, 0, 0, 0, damp);
    vFE->SetPlane(0, -1, 0, 0, 0, 0.2);
    vFE->SetFillMesh(false);
    const bool plotDrift2 = true;
    if (plotDrift2) {
        TCanvas* c2 = new TCanvas();
        vFE->SetCanvas(c2);
        //for (unsigned int i = 0; i < nMaterials; ++i) vFE->SetColor(i, kGray+i);
        //vFE->SetColor(0, kWhite);
        //vFE->SetColor(0, kYellow + 3);
        //vFE->SetColor(1, kGray);
        //vFE->SetColor(2, kYellow + 3);
        vFE->EnableAxes();
        vFE->SetXaxisTitle("x (cm)");
        vFE->SetYaxisTitle("z (cm)");
        std::cout << "Plotting..." << std::endl;
        vFE->Plot();
        c2->SaveAs(Form("Figures/avalanche2d_%s_model%d.pdf", gasName.c_str(), modelNum));
    }
    
    time_t t1 = time(NULL);
    
    //std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
    
    return 0;
}


