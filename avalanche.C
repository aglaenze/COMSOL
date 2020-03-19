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
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 8;
    //____________________
    
    time_t t0 = time(NULL);
    
    // Make a gas medium.
    MediumMagboltz* gas = InitiateGas(gasName);
    // Load field map
    ComponentComsol* fm;
    
    TString fOutputName;
    int hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvGemDown = 0, hvGemUp = 0, hvDrift = 0;
    if (modelNum == 1) {
        if (argc < 3) {
            std::cout << "Please enter command like this: ./avalanche $hvMesh $hvDrift " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDrift = atoi(argv[2]);
        fm = InitiateField(modelNum, hvMesh, hvDrift, gas);
        fOutputName = Form("Figures/avalanche2d-%s-model%d-%d-%d.pdf", gasName.c_str(), modelNum, hvMesh, hvDrift);
    }
    else if (modelNum >= 2 && modelNum < 5) {
        if (argc < 4) {
            std::cout << "Please enter command like this: ./avalanche $hvDmDown $hvDmUp $hvDrift " << std::endl;
            return 0;
        }
        hvDmDown = atoi(argv[1]);
        hvDmUp = atoi(argv[2]);
        hvDrift = atoi(argv[3]);
        fm = InitiateField(modelNum, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("Figures/avalanche2d-%s-model%d-%d-%d-%d.pdf", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift);
    }
    else if (modelNum >= 5 && modelNum < 8) {
        if (argc < 5) {
            std::cout << "Please enter command like this: ./avalanche $hvMesh $hvDmDown $hvDmUp $hvDrift " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDmDown = atoi(argv[2]);
        hvDmUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        fm = InitiateField(modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("Figures/avalanche2d-%s-model%d-%d-%d-%d-%d.pdf", gasName.c_str(), modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift);
    }
    else if (modelNum >= 8 && modelNum < 10) {
        if (argc != 5) {
            std::cout << "Please enter HVmesh like this: ./avalanche $hvMesh $hvGemDown $hvGemUp $hvDrift " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvGemDown = atoi(argv[2]);
        hvGemUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        fm = InitiateField(modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift, gas);
        fOutputName = Form("Figures/avalanche2d-%s-model%d-%d-%d-%d-%d.pdf", gasName.c_str(), modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift);
    }
    else {std::cout << "Wrong model number" << std::endl; return 0;}
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
    
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
    drift->SetDistanceSteps(2.e-4);
    
    aval->EnablePlotting(driftView);
    aval->EnablePlotting(driftView2);
    drift->EnablePlotting(driftView);
    drift->EnablePlotting(driftView2);
    
    
    const int nEvents = 1;
    
    int j = 0;
    for (unsigned int i = 0; i < nEvents; ++i) {
        j++;
        if (j==50) break;
        //std::cout << "hello " << i << "\n\n" << std::endl;
        // Initial coordinates of the electron.
        double x0 = RndmUniform() * pitch;
        double y0 = RndmUniform() * depth;
        double z0 = damp + 2*radius + (ddrift-damp-2*radius)*RndmUniform();
        //double z0 = 0.3;
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
            //drift->DriftIon(xe1, ye1, ze1, te1);
            //drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
        }
    }
    
    //return 0;
    
    const bool plotDrift = false;
    if (plotDrift) {
        TCanvas* c1 = new TCanvas();
        driftView->SetCanvas(c1);
        driftView->Plot();
        c1->SaveAs(Form("Figures/avalanche-%s.pdf", gasName.c_str()));
    }
  
        
    // Set up the object for FE mesh visualization.
    ViewFEMesh* vFE = new ViewFEMesh();
    //vFE->SetArea(-5*pitch, -5*pitch, 0,  5*pitch, 5*pitch, damp*2);
    //vFE->SetArea(-5*pitch, -5*pitch, damp-5*pitch,  5*pitch, 5*pitch, damp+5*pitch);
    if (modelNum > 1 && modelNum < 5) {vFE->SetArea(-25*pitch, -25*pitch, 0,  25*pitch, 25*pitch, 50*pitch);}
    else if (modelNum >= 5 && modelNum < 8) {vFE->SetArea(-damp/2, -damp/2, 0,  damp/2+2*pitch, damp/2+2*pitch, damp+2*pitch);}
    //vFE->SetArea(-10*pitch, -10*pitch, damp-15*pitch,  10*pitch, 10*pitch, damp+5*pitch);
    vFE->SetArea(-(damp+5*pitch)/2., -(damp+5*pitch)/2., 0,  (damp+5*pitch)/2., (damp+5*pitch)/2., damp+5*pitch);
    vFE->SetComponent(fm);
    vFE->SetViewDrift(driftView2);
    //vFE->SetPlane(0, -1, 0, 0, 0, damp);
    //vFE->SetPlane(0, -1, 0, 0, 0, 0);
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
        c2->SaveAs(fOutputName);
    }
    
    time_t t1 = time(NULL);
    
    //std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
    
    return 0;
}


