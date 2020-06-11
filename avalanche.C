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
	int modelNum = 0;
	std::string gasName = "";
	if(!LoadVariables(modelNum, gasName)) {std::cout << "variables not loaded" << std::endl; return 0;}
    const bool plotDrift3D = false;
	const bool plotDrift2D = true;
    //____________________
    
	time_t t0 = time(NULL);
	if (modelNum < 1 || modelNum > GetMaxModelNum()) {std::cout << "Wrong model number" << std::endl; return 0;}
	
	TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();
	
	int electrodeNum = 0;
	electrodeNum = GetElectrodeNum(modelNum);
	if (electrodeNum == 0) {std::cout << "Warning! Number of electrodes = 0" << std::endl; return 0;}
	
	TString errorMessage = "Please enter HVmesh like this: ./avalanche";
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
	
	TString fOutputName = Form("Figures/model%d/avalanche2d-%s", modelNum,  gasName.c_str());
	for (int k = 0; k< electrodeNum-1; k++) fOutputName += Form("-%d", hvList[k]);
	fOutputName += ".pdf";
	
	//Load geometry parameters
	double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	int periodicityNum = 0;
	LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
    
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
    
	// To look at the avalanche
	ViewDrift* driftView = new ViewDrift();
	driftView->SetPlane(0, -1, 0, 0, 0, 0);
    aval->EnablePlotting(driftView);
    if (plotDrift3D) drift->EnablePlotting(driftView);
    
    const int nEvents = 1;
    
    int j = 0;
    for (unsigned int i = 0; i < nEvents; ++i) {
        j++;
        if (j==50) break;
        //std::cout << "hello " << i << "\n\n" << std::endl;
        // Initial coordinates of the electron.
        double x0 = 0;
        double y0 = 0;
        //double z0 = damp + 2*radius + (ddrift-damp-2*radius)*RndmUniform();
        double z0 = 0.03;
        double t0 = 0;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        int ne2 = 0, ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        std::cout << "\nAvalanche size = " << ne2 << std::endl;
		std::cout << "\nIon avalanche size = " << ni << std::endl;
        if (ne2 < 100) {i--; continue;}
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
		driftView->SetArea(-10*pitch, -10*pitch, 0, 10*pitch, 10*pitch, damp*5);
        driftView->SetCanvas(c1);
        driftView->Plot(false, true);
        c1->SaveAs(Form("Figures/avalanche-%s.pdf", gasName.c_str()));
    }
  

    // Set up the object for FE mesh visualization.
    ViewFEMesh* vFE = new ViewFEMesh();
    //vFE->SetArea(-5*pitch, -5*pitch, 0,  5*pitch, 5*pitch, damp*2);
    //vFE->SetArea(-5*pitch, -5*pitch, damp-5*pitch,  5*pitch, 5*pitch, damp+5*pitch);
    if ((modelNum > 1 && modelNum < 5) || modelNum == 14) {vFE->SetArea(-25*pitch, -25*pitch, 0,  25*pitch, 25*pitch, 50*pitch);}
    else if (modelNum >= 5 && modelNum < 8) {vFE->SetArea(-damp/2, -damp/2, 0,  damp/2+2*pitch, damp/2+2*pitch, damp+2*pitch);}
	else {vFE->SetArea(-(damp+5*pitch)/2., -(damp+5*pitch)/2., 0,  (damp+5*pitch)/2., (damp+5*pitch)/2., damp+5*pitch);}
    vFE->SetComponent(fm);
	driftView->SetArea(-5*pitch, -5*pitch, 0, 5*pitch, 5*pitch, damp*2);
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


