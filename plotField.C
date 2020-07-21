#include <iostream>
#include <fstream>
#include <ctime>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "initiate.C"
#include "_Utils.C"
#include "_Geometry.C"

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
    
	// variables, to change in the file input.txt
	int modelNum = 0;
	std::string gasName = "";
	if(!LoadVariables(modelNum, gasName)) {std::cout << "variables not loaded" << std::endl; return 0;}
    //____________________
	//time_t t0 = time(NULL);
	if (modelNum < 1 || modelNum > GetMaxModelNum()) {std::cout << "Wrong model number" << std::endl; return 0;}
	
	TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();
	
	int electrodeNum = 0;
	electrodeNum = GetElectrodeNum(modelNum);
	if (electrodeNum == 0) {std::cout << "Warning! Number of electrodes = 0" << std::endl; return 0;}
	
	TString errorMessage = "Please enter HVmesh like this: ./plotField";
	for (int k = 0; k< electrodeNum; k++) errorMessage += Form(" $hv%d", k+1);
	if (argc != electrodeNum) {
		std::cout << errorMessage << std::endl;
		return 0;
	}
	std::vector<int> hvList = {};
	for (int k = 1; k < electrodeNum; k++) hvList.push_back(atoi(argv[k]) );
	
	// Make a gas medium.
	MediumMagboltz* gas = InitiateGas(gasName);
	ComponentComsol* fm = InitiateField(modelNum, hvList, gas);
	if (!fm || fm->GetMedium(0,0,0) == nullptr) {
		std::cout << "Component COMSOL was not initialized, please fix this" << std::endl;
		return 0;
	}
	
	//Load geometry parameters
	double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    ViewField* vf = new ViewField();
    vf->SetComponent(fm);
    //if (modelNum == 4) {vf->SetPlane(1, -1, 0, 0, 0, 0); vf->Rotate(TMath::Pi()*1.5);}
    //else if (modelNum == 5) vf->SetPlane(0, -1, 0, 0, 0.5*pitch, 0);
    //else vf->SetPlane(0, -1, 0, 0, 0, 0);
    vf->SetPlane(0, -1, 0, 0, 0, 0);
    if (modelNum>=8 && modelNum<14) vf->SetPlane(0, -1, 0, pitch/2, pitch/2, 0);
    //vf->SetVoltageRange(-550, 0);
    
    const bool plotField = false;
    if (plotField) {
        //vf->SetVoltageRange(-600., 0.);
        //vf->SetArea(0, 0, width, 3*damp);
        //vf->SetArea(0, 0, width, 3*damp+0.2);
        vf->SetArea(0, 0, width, ddrift);
        TCanvas* c1 = new TCanvas("c1", "Potential view", 1200, 600);
		TCanvas* c2 = new TCanvas("c2", "Field view", 1200, 600);
        //TCanvas* cf = new TCanvas("cf", "Potential view", 600, 600);
        vf->SetCanvas(c1);
		vf->PlotContour("v");
        c1->SaveAs(Form("Figures/model%d/potential.pdf", modelNum));
		vf->SetCanvas(c2);
        vf->PlotContour("e");
        c2->SaveAs(Form("Figures/model%d/field.pdf", modelNum));

    }
    
    vf->SetVoltageRange(-hvList[1], -hvList[2]);
    const bool zoom = true;
    if (zoom) {
        vf->SetArea(0, damp-4*pitch, 5*pitch, damp+pitch);
        vf->SetArea(0, damp-pitch, 2*pitch, damp+pitch);
        //vf->SetNumberOfContours(40);
        //vf->SetNumberOfSamples2d(40, 40);
		TCanvas* c1 = new TCanvas("c11", "c11", 600, 600);
        TCanvas* c2 = new TCanvas("c22", "c22", 600, 600);
		TCanvas* c3 = new TCanvas("c3", "Field lines", 1200, 600);
        //TCanvas* c2 = new TCanvas("c2", "c2", 1000*4*pitch, 1000*damp);
        vf->SetCanvas(c1);
        if (modelNum==1) vf->SetVoltageRange(-hvList[0]*1.1, -hvList[0]*0.78);
        //else if (modelNum>=8 && modelNum<14) vf->SetVoltageRange(-hvList[1]*1.05, -hvList[1]/1.05);
        vf->PlotContour("v");
        c1->SaveAs(Form("Figures/model%d/potentialZoom.pdf", modelNum));
		vf->SetCanvas(c2);
        vf->PlotContour("e");
        c2->SaveAs(Form("Figures/model%d/fieldZoom.pdf", modelNum));
		vf->SetCanvas(c3);
		//vf->SetNumberOfContours(2);
		//vf->PlotProfile(0, 0, 0, width/2, width/2, damp*2, "v");
		vf->PlotContour("ey");
		c3->SaveAs(Form("Figures/model%d/ey.pdf", modelNum));
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
        c3->SaveAs(Form("Figures/model%d/Mesh.pdf", modelNum));
    }
    

    //app.Run(kTRUE);
    
}

