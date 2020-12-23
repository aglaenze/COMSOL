#include <iostream>
#include <fstream>
#include <ctime>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "Include/Initiate.C"
#include "Include/Utils.C"
#include "Include/Geometry.C"

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
using namespace std;

int main(int argc, char * argv[]) {
	
	// variables, to change in the file input.txt
	int modelNum = 0;
	std::string gasName = "";
	if(!LoadVariables(modelNum, gasName)) {std::cout << "variables not loaded" << std::endl; return 0;}
	bool remote = false;
	//____________________
	//time_t t0 = time(NULL);
	if (modelNum < 1 || modelNum > GetMaxModelNum(remote)) {std::cout << "Wrong model number" << std::endl; return 0;}
	
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
	vector<int> hvList = {};
	for (int k = 1; k < electrodeNum; k++) hvList.push_back(atoi(argv[k]) );
	
	TString suffix = "";
	for (int k = 0; k< electrodeNum-1; k++) {suffix += Form("-%d", hvList[k]);}
	suffix += ".pdf";
	
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
	
	// Make a sensor.
	Sensor sensor;
	sensor.AddComponent(fm);
	sensor.SetArea(-width/2, -depth/2, 0, width/2, depth/2, ddrift);
	
	ViewField* vf = new ViewField();
	//vf->SetComponent(fm);
	vf->SetSensor(&sensor);
	vf->SetArea(0, 0, width, ddrift);
	double yPlane = pitch/4;
	//double yPlane = 0;
	vf->SetPlane(0, -1, 0, 0, yPlane, 0);
	//vf->SetPlaneXZ();
	/*
	vf->SetPlane(1, -1, 0, pitch/4, pitch/4, 0);
	vf->Rotate(TMath::Pi()*1.5);
	 */
	if (electrodeNum>3) vf->SetVoltageRange(-hvList[0], -hvList[electrodeNum-2]);
	
	const bool plotField = true;
	if (plotField) {
		//vf->SetNumberOfContours(10);
		//vf->SetNumberOfSamples2d(40, 40);
		TCanvas* c1 = new TCanvas("c1", "Potential view", 1200, 600);
		TCanvas* c2 = new TCanvas("c2", "Field view", 1200, 600);
		//TCanvas* cf = new TCanvas("cf", "Potential view", 600, 600);
		vf->SetCanvas(c1);
		vf->Plot("v", "CONT4Z");
		DrawElectrodes(modelNum, 0, ddrift);
		c1->SaveAs(Form("Figures/model%d/potential", modelNum)+suffix);
		vf->SetCanvas(c2);
		vf->Plot("e", "CONT4Z");
		gPad->Update();
		DrawElectrodes(modelNum, 0, ddrift);
		c2->SaveAs(Form("Figures/model%d/field", modelNum)+suffix);
	}
	
	const bool zoom = true;
	if (zoom) {
		/*
		double size = 4*pitch;	// size of a side of the window to plot
		const double xmin = 0;
		double zmin = damp-size*0.7;
		if (zmin < 0) {zmin = 0;}
		if (zmin+size > ddrift) {size = ddrift-zmin;}
		const double xmax =  xmin+size;
		const double zmax =  damp+size*0.3;
		 */
	
		double zmin = damp*0.9, zmax = damp*1.05;
		double xmin = -(zmax-zmin)/2, xmax = (zmax-zmin)/2;
		vf->SetArea(xmin, zmin, xmax, zmax);
		//vf->SetArea(0, damp-pitch, 2*pitch, damp+pitch);
		TCanvas* c1 = new TCanvas("c11", "c11", 600, 600);
		TCanvas* c2 = new TCanvas("c22", "c22", 600, 600);
		//TCanvas* c2 = new TCanvas("c2", "c2", 1000*4*pitch, 1000*damp);
		vf->SetCanvas(c1);
		if (modelNum==1) vf->SetVoltageRange(-hvList[0]*1.1, -hvList[0]*0.78);
		//else if (modelNum>=8 && modelNum<14) vf->SetVoltageRange(-hvList[1]*1.05, -hvList[1]/1.05);
		vf->Plot("v", "CONT4Z");
		DrawElectrodes(modelNum, zmin, zmax);
		c1->SaveAs(Form("Figures/model%d/potentialZoom", modelNum)+suffix);
		vf->SetCanvas(c2);
		vf->Plot("e", "CONT4Z");
		//DrawElectrodes(modelNum, zmin, zmax);
		c2->SaveAs(Form("Figures/model%d/fieldZoom", modelNum)+suffix);
		
		// Field lines
		//TCanvas* c3 = new TCanvas("c3", "Field lines", 600, 600);
		zmin = 0;
		zmax = ddrift;
		if ((modelNum > 1 && modelNum < 5) || modelNum == 14) {zmax = 50*pitch;}
		else if (modelNum >= 5 && modelNum < 8) {zmax = damp+2*pitch;}
		else {zmax = damp+5*pitch;}
		zmin = damp*0.9; zmax = damp*1.05;
		xmin = -(zmax-zmin)/2., xmax = (zmax-zmin)/2.;
		vf->SetArea(xmin, zmin, xmax, zmax);
		
		TCanvas* c3 = new TCanvas();
		vf->SetCanvas(c3);
		vector<double> xf;
		vector<double> yf;
		vector<double> zf;
		int nPitch = int((xmax-xmin)/pitch);
		//cout << "nPitch = " << nPitch << endl;
		int nLines = nPitch*20-1;
		double xLineMin = -nPitch*pitch/2;
		double xLineMax = nPitch*pitch/2;
		vf->EqualFluxIntervals(xLineMin, yPlane, zmax*0.99, xLineMax, yPlane, zmax*0.99, xf, yf, zf, nLines);
		
		//vf->EqualFluxIntervals(xmin, -pitch, 0.99 * zmax, xmax, pitch, 0.99 * zmax, xf, yf, zf, 200);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetRightMargin(0.15);
		//vf->Plot("v", "CONT4Z");	// "CONT1Z"
		//vf->Plot("v", "CONT1");
		vf->PlotFieldLines(xf, yf, zf, true, true);	// last one should be false in you want to plot something else before (pltaxis = false)
		//DrawElectrodes(modelNum, zmin, zmax);
		c3->SaveAs(Form("Figures/model%d/fieldlines", modelNum)+suffix);
		//vf->SetNumberOfContours(2);
		//vf->PlotProfile(0, 0, 0, width/2, width/2, damp*2, "v");
		/*
		TCanvas* c4 = new TCanvas("c4", "Transverse field", 600, 600);
		vf->SetCanvas(c4);
		//vf->PlotContour("ey");
		vf->Plot("ey", "CONT4Z");
		DrawElectrodes(modelNum, zmin, zmax);
		c4->SaveAs(Form("Figures/model%d/ey", modelNum)+suffix);
		 */
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
		c3->SaveAs(Form("Figures/model%d/Mesh", modelNum)+suffix);
	}
	
	
	//app.Run(kTRUE);
	
}

