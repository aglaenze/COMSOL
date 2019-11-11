#include <iostream>
#include <fstream>
#include <ctime>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>

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
    const int modelNum = 1;
    //____________________
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    if (argc < 2) {
        std::cout << "Please enter a HVmesh like this: ./plotField $hvMesh " << std::endl;
        return 0;
    }
    const int hvMesh = atoi(argv[1]);

      // Make a gas medium.
      MediumMagboltz* gas = InitiateGas(gasName);
      // Load field map
      ComponentComsol* fm = InitiateField(modelNum, hvMesh, gas);

    
    const bool plotField = true;
    if (plotField) {
        ViewField* fieldView = new ViewField();
        fieldView->SetComponent(fm);
        //fieldView->SetVoltageRange(-600., 0.);
        fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
        fieldView->SetArea(0, 0, width, 3*damp);
        //fieldView->SetVoltageRange(-hvMesh*1.2, 0);
        TCanvas* c1 = new TCanvas("cf", "Potential view", 1200, 600);
        //TCanvas* cf = new TCanvas("cf", "Potential view", 600, 600);
        fieldView->SetCanvas(c1);
        fieldView->PlotContour();
        c1->SaveAs(Form("Figures/potential_model%d.pdf", modelNum));
    }
    
    const bool zoom = true;
    if (zoom) {
        ViewField* vf = new ViewField();
        vf->SetComponent(fm);
        //vf->SetPlane(0, -1, 0, 0, 0.5*pitch, 0);
        vf->SetPlane(0, -1, 0, 0, 0, 0);
        vf->SetArea(pitch, damp-pitch, 3*pitch, damp+pitch);
        vf->SetNumberOfContours(70);
        vf->SetNumberOfSamples2d(40, 40);
        TCanvas* c2 = new TCanvas();
        c2->SetLeftMargin(0.2);
        vf->SetCanvas(c2);
        vf->SetVoltageRange(-hvMesh*1.1, -hvMesh*0.78);
        //vf->PlotContour("v");
        vf->PlotContour("v");
        //c2->SaveAs("Figures/potentialZoom.pdf");
        c2->SaveAs(Form("Figures/fieldZoom_model%d.pdf", modelNum));
    }
    

    const bool plotMesh = true;
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

