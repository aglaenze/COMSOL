#include <iostream>
#include <fstream>
#include <ctime>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "parameters.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/ComponentBase.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {
    
    // variables
    std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 6;
    //____________________
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    if (argc < 2) {
        std::cout << "Please enter a HVmesh like this: ./plotField $hvMesh " << std::endl;
        return 0;
    }
    const int hvMesh = atoi(argv[1]);

    //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
    
      // Make a gas medium.
      MediumMagboltz* gas = InitiateGas(gasName);
      // Load field map
      ComponentComsol* fm = InitiateField(modelNum, hvMesh, gas);
    //return 0;
    
    ViewField* vf = new ViewField();
    vf->SetComponent(fm);
    if (modelNum == 4) {vf->SetPlane(1, -1, 0, 0, 0, 0); vf->Rotate(TMath::Pi()*1.5);}
    else if (modelNum == 5) vf->SetPlane(0, -1, 0, 0, 0.5*pitch, 0);
    else vf->SetPlane(0, -1, 0, 0, 0, 0);
    //vf->SetVoltageRange(-hvMesh*1., 0);
    vf->SetVoltageRange(-550, 0);
    
    const bool plotField = true;
    if (plotField) {
        //vf->SetVoltageRange(-600., 0.);
        //vf->SetArea(0, 0, width, 3*damp);
        vf->SetArea(0, 0, width, 3*damp+0.2);
        TCanvas* c1 = new TCanvas("cf", "Potential view", 1200, 600);
        //TCanvas* cf = new TCanvas("cf", "Potential view", 600, 600);
        vf->SetCanvas(c1);
        vf->PlotContour();
        c1->SaveAs(Form("Figures/potential_model%d.pdf", modelNum));
        vf->PlotContour("e");
        c1->SaveAs(Form("Figures/field_model%d.pdf", modelNum));
    }
    
    const bool zoom = true;
    if (zoom) {
        //vf->SetArea(pitch, damp-pitch, 3*pitch, damp+pitch);
        //vf->SetArea(pitch, 0, 5*pitch, damp+pitch);
        vf->SetArea(pitch, 0.2120-pitch, 5*pitch, 0.2120+4*pitch);
        vf->SetNumberOfContours(70);
        vf->SetNumberOfSamples2d(40, 40);
        TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
        //TCanvas* c2 = new TCanvas("c2", "c2", 1000*4*pitch, 1000*damp);
        c2->SetLeftMargin(0.1);
        vf->SetCanvas(c2);
        //vf->SetVoltageRange(-hvMesh*1.1, -hvMesh*0.78);
        vf->PlotContour("e");
        c2->SaveAs(Form("Figures/fieldZoom_model%d.pdf", modelNum));
        vf->PlotContour("v");
        c2->SaveAs(Form("Figures/potentialZoom_model%d.pdf", modelNum));
    }
    

    const bool plotMesh = false;
    if (plotMesh) {
        TCanvas* c3 = new TCanvas();
        ViewFEMesh* meshView = new ViewFEMesh();
        // Set the component.
        meshView->SetComponent(fm);
        // Set the viewing plane.
        //meshView->SetPlane(0., -1., 0., 2.5*pitch, 2.5*pitch, damp);
        meshView->SetPlane(0., -1., 0., 2*pitch, 2*pitch, damp);
        meshView->SetColor(0, 1);     // matid, colorid
        //meshView->SetFillColor(0, 0);
        meshView->SetColor(1, 4);     // matid, colorid
        //meshView->SetFillColor(1, 4);
        //meshView->SetFillMesh(true);
        meshView->SetCanvas(c3);
        meshView->SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
        meshView->EnableAxes();
        meshView->Plot();
        c3->SaveAs(Form("Figures/Mesh_model%d.pdf", modelNum));
    }
    

    //app.Run(kTRUE);
    
}

