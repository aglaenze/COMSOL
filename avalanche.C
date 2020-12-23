#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "Include/Utils.C"
#include "Include/Geometry.C"
#include "Include/Initiate.C"

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
using namespace std;

int main(int argc, char * argv[]) {
	
	//______________________
	// variables
	int modelNum = 0;
	string gasName = "";
	bool remote = false;
	bool plotDrift2D = 0, plotDrift3D = 0;
	if(!LoadVariables(modelNum, gasName, plotDrift2D, plotDrift3D)) {cout << "variables not loaded" << endl; return 0;}
	//____________________
	
	time_t t0 = time(NULL);
	if (modelNum < 1 || modelNum > GetMaxModelNum(remote)) {cout << "Wrong model number" << endl; return 0;}

	TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();
	
	int electrodeNum = 0;
	electrodeNum = GetElectrodeNum(modelNum);
	if (electrodeNum == 0) {cout << "Warning! Number of electrodes = 0" << endl; return 0;}
	
	TString errorMessage = "Please enter HVmesh like this: ./avalanche";
	for (int k = 0; k< electrodeNum; k++) errorMessage += Form(" $hv%d", k+1);
	if (argc != electrodeNum) {
		cout << errorMessage << endl;
		return 0;
	}
	vector<int> hvList = {};
	for (int k = 1; k < electrodeNum; k++) hvList.push_back(atoi(argv[k]) );
	
	// Make a gas medium.
	MediumMagboltz* gas = InitiateGas(gasName);
	ComponentComsol* fm = InitiateField(modelNum, hvList, gas);
	if (!fm || fm->GetMedium(0,0,0) == nullptr) {
		cout << "Component COMSOL was not initialized, please fix this" << endl;
		return 0;
	}
	
	TString fOutputName2d = Form("Figures/model%d/avalanche2d-%s", modelNum,  gasName.c_str());
	TString fOutputName2dZoom = Form("Figures/model%d/avalanche2d-zoom-%s", modelNum,  gasName.c_str());
	TString fOutputNameIons2d = Form("Figures/model%d/avalanche2d-ions-%s", modelNum,  gasName.c_str());
	TString fOutputNameIons2dZoom = Form("Figures/model%d/avalanche2d-ions-zoom-%s", modelNum,  gasName.c_str());
	TString fOutputName3d = Form("Figures/model%d/avalanche3d-%s", modelNum,  gasName.c_str());
	for (int k = 0; k< electrodeNum-1; k++) {
		fOutputName2d += Form("-%d", hvList[k]);
		fOutputName2dZoom += Form("-%d", hvList[k]);
		fOutputNameIons2d += Form("-%d", hvList[k]);
		fOutputNameIons2dZoom += Form("-%d", hvList[k]);
		fOutputName3d += Form("-%d", hvList[k]);
	}
	fOutputName2d += ".pdf";
	fOutputName2dZoom += ".pdf";
	fOutputNameIons2d += ".pdf";
	fOutputNameIons2dZoom += ".pdf";
	fOutputName3d += ".pdf";
	
	//Load geometry parameters
	double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
	
	// Make a sensor.
	Sensor sensor;
	sensor.AddComponent(fm);
	sensor.SetArea(-width/2, -depth/2, 0, width/2, depth/2, ddrift);
	//sensor.SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
	
	
	// Create an avalanche object
	AvalancheMicroscopic* aval = new AvalancheMicroscopic();
	aval->SetSensor(&sensor);
	AvalancheMC* drift = new AvalancheMC();
	drift->SetSensor(&sensor);
	drift->SetDistanceSteps(2.e-4);
	
	// To look at the avalanche of electrons (if 3D, ions too)
	ViewDrift* driftView = new ViewDrift();
	driftView->SetPlane(0, -1, 0, 0, 0, 0);
	aval->EnablePlotting(driftView);
	if (plotDrift3D) drift->EnablePlotting(driftView);
	
	// To look at the avalanche of ions
	ViewDrift* driftViewIons = new ViewDrift();
	driftViewIons->SetPlane(0, -1, 0, 0, 0, 0);
	drift->EnablePlotting(driftViewIons);
	
	const int nEvents = 1;
	
	int j = 0;
	for (unsigned int i = 0; i < nEvents; ++i) {
		j++;
		if (j==50) break;
		//cout << "hello " << i << "\n\n" << endl;
		// Initial coordinates of the electron.
		double x0 = 0;
		double y0 = 0;
		//double z0 = damp + 2*radius + (ddrift-damp-2*radius)*RndmUniform();
		double z0 = damp+ (ddrift-damp)/2;
		double t0 = 0;
		double e = 0;
		aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
		int ne2 = 0, ni = 0;
		aval->GetAvalancheSize(ne2, ni);
		cout << "\nAvalanche size = " << ne2 << endl;
		cout << "\nIon avalanche size = " << ni << endl;
		if (ne2 < 100) {i--; driftView->Clear(); continue;}
		const int np = aval->GetNumberOfElectronEndpoints();
		double xe1, ye1, ze1, te1, e1;
		double xe2, ye2, ze2, te2, e2;
		double xi1, yi1, zi1, ti1;
		double xi2, yi2, zi2, ti2;
		int status;
		for (int j = np; j--;) {
			aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
			if (plotDrift3D) {
				drift->DriftIon(xe1, ye1, ze1, te1);
				drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
			}
		}
		
		
	}
	
	//return 0;
	
	if (plotDrift3D) {
		TCanvas* c1 = new TCanvas();
		driftView->SetArea(-5*pitch, -5*pitch, 0, 5*pitch, 5*pitch, damp*2);
		driftView->SetCanvas(c1);
		driftView->Plot(false, true);
		c1->SaveAs(fOutputName3d);
	}
	
	
	// Set up the object for FE mesh visualization.
	ViewFEMesh* vFE = new ViewFEMesh();
	//vFE->SetArea(-5*pitch, -5*pitch, 0,  5*pitch, 5*pitch, damp*2);
	//vFE->SetArea(-5*pitch, -5*pitch, damp-5*pitch,  5*pitch, 5*pitch, damp+5*pitch);
	double zmin = 0;
	double zmax = ddrift;
	if ((modelNum > 1 && modelNum < 5) || modelNum == 14) {zmax = 50*pitch;}
	else if (modelNum >= 5 && modelNum < 8) {zmax = damp+2*pitch;}
	else {zmax = damp+5*pitch;}
	
	vFE->SetComponent(fm);
	driftView->SetArea(-5*pitch, -5*pitch, 0, 5*pitch, 5*pitch, damp*2);
	//vFE->SetArea(-5*pitch, -5*pitch, damp-8*pitch, 5*pitch, 5*pitch, damp+2*pitch);
	vFE->SetViewDrift(driftView);
	vFE->SetPlane(0, -1, 0, 0, 0, 0);
	vFE->SetFillMesh(true);
	if (plotDrift2D) {
		TCanvas* c2 = new TCanvas();
		vFE->SetCanvas(c2);
		/*
		 //for (unsigned int i = 0; i < nMaterials; ++i) vFE->SetColor(i, kGray+i);
		 //vFE->SetColor(0, kWhite);
		 vFE->SetColor(0, kYellow + 3);
		 vFE->SetColor(1, kGray);
		 vFE->SetColor(2, kGray);
		 */
		vFE->SetArea(-(zmax-zmin)/2., -(zmax-zmin)/2., zmin,  (zmax-zmin)/2., (zmax-zmin)/2., zmax);
		vFE->EnableAxes();
		vFE->SetXaxisTitle("x (cm)");
		vFE->SetYaxisTitle("z (cm)");
		cout << "Plotting..." << endl;
		vFE->Plot();
		DrawElectrodes(modelNum, zmin, zmax);
		c2->SaveAs(fOutputName2d);
		
		zmin = damp*0.9; zmax = damp*1.05;
		vFE->SetArea(-(zmax-zmin)/2., -(zmax-zmin)/2., zmin,  (zmax-zmin)/2., (zmax-zmin)/2., zmax);
		TCanvas* c3 = new TCanvas();
		vFE->SetCanvas(c3);
		cout << "Plotting..." << endl;
		vFE->Plot();
		DrawElectrodes(modelNum, zmin, zmax);
		c3->SaveAs(fOutputName2dZoom);
		
		// Same with ions
		vFE->SetViewDrift(driftViewIons);
		vFE->SetArea(-(zmax-zmin)/2., -(zmax-zmin)/2., zmin,  (zmax-zmin)/2., (zmax-zmin)/2., zmax);
		TCanvas* c4 = new TCanvas();
		vFE->SetCanvas(c4);
		cout << "Plotting..." << endl;
		vFE->Plot();
		DrawElectrodes(modelNum, zmin, zmax);
		c4->SaveAs(fOutputNameIons2d);
		
		zmin = damp*0.9; zmax = damp*1.05;
		vFE->SetArea(-(zmax-zmin)/2., -(zmax-zmin)/2., zmin,  (zmax-zmin)/2., (zmax-zmin)/2., zmax);
		TCanvas* c5 = new TCanvas();
		vFE->SetCanvas(c5);
		cout << "Plotting..." << endl;
		vFE->Plot();
		DrawElectrodes(modelNum, zmin, zmax);
		c5->SaveAs(fOutputNameIons2dZoom);
	}
	
	
	time_t t1 = time(NULL);
	
	//cout << "\n" << nEvents << " events simulated" << endl;
	PrintTime(t0, t1);
	
	if (plotDrift3D) app.Run(true);
	
	return 0;
}


