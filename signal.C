#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "_Utils.C"
#include "_Geometry.C"
#include "initiate.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {
	
	bool testMode = false;
	bool keepSignal = false;
	//______________________
	// variables, to change in the file input.txt
	int modelNum = 0;
	string gasName = "";
	bool computeIBF = true;
	bool useFeSource = true;
	int nEvents = 0;  // number of avalanches to simulate
	if(!LoadVariables(modelNum, gasName, nEvents, computeIBF, useFeSource)) {cout << "variables not loaded" << endl; return 0;}
	//____________________
	
	
	time_t t0 = time(NULL);
	if (modelNum < 1 || modelNum > GetMaxModelNum()) {cout << "Wrong model number" << endl; return 0;}
	
	TApplication app("app", &argc, argv);
	plottingEngine.SetDefaultStyle();
	
	int electrodeNum = 0;
	electrodeNum = GetElectrodeNum(modelNum);
	if (electrodeNum == 0) {cout << "Warning! Number of electrodes = 0" << endl; return 0;}
	
	TString errorMessage = "Please enter HVmesh like this: ./signal";
	for (int k = 0; k< electrodeNum; k++) errorMessage += Form(" $hv%d", k+1);
	if (testMode && argc != electrodeNum) {
		errorMessage += " (test mode)";
		cout << errorMessage << endl;
		return 0;
	}
	if (!testMode && argc != electrodeNum+1) {
		errorMessage += " $saveNum";
		cout << errorMessage << endl;
		return 0;
	}
	vector<int> hvList = {};
	for (int k = 1; k < electrodeNum; k++) hvList.push_back(atoi(argv[k]) );
	int saveNum = atoi(argv[electrodeNum]);
	
	// Make a gas medium.
	MediumMagboltz* gas = InitiateGas(gasName);
	ComponentComsol* fm = InitiateField(modelNum, hvList, gas);
	if (!fm || fm->GetMedium(0,0,0) == nullptr) {
		cout << "Component COMSOL was not initialized, please fix this" << endl;
		return 0;
	}
	
	string type = "signal";
	if (useFeSource) type = "fe-signal";
	if (!computeIBF) type += "-noibf";
	TString fOutputName = Form("rootFiles/%s/model%d/%s", gasName.c_str(), modelNum, type.c_str());
	for (int k = 0; k< electrodeNum-1; k++) fOutputName += Form("-%d", hvList[k]);
	if (testMode) fOutputName += "-test.root";
	else fOutputName += Form("-%d.root", saveNum);
	
	//Load geometry parameters
	double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
	
	
	// Make a sensor.
	Sensor* sensor = new Sensor();
	sensor->AddComponent(fm);
	sensor->SetArea(0, 0, 0, width, depth, ddrift);
	//sensor->SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
	for (int k = 0; k < electrodeNum; k++) sensor->AddElectrode(fm, Form("V%d", k+2));
	//return 0;
	
	// Use Heed for simulating the photon absorption.
	TrackHeed track;
	if (useFeSource) {
		track.SetSensor(sensor);
		track.EnableElectricField();
	}
	
	// Create ROOT histograms of the signal and a file in which to store them.
	TFile* f = new TFile(fOutputName, "RECREATE");
	//TFile* f = new TFile("rootFiles/test.root", "RECREATE");
	
	TTree *tAvalanche = new TTree("tAvalanche","Gain");
	int nWinners = 0, neAval = 0;
	vector<int> neAvalVec = {}, nWinnersVec = {};
	tAvalanche->Branch("amplificationElectrons", &nWinnersVec);
	tAvalanche->Branch("avalancheSize", &neAvalVec);
	int ni = 0, ionBackNum = 0;
	vector<float> electronStartPoints = {}, electronEndPoints = {};
	vector<float> ionStartPoints = {}, ionEndPoints = {}, ionEndPointsX = {}, ionEndPointsY = {};
	if (computeIBF) {
		tAvalanche->Branch("electronStartPoints", &electronStartPoints);
		tAvalanche->Branch("electronEndPoints", &electronEndPoints);
		tAvalanche->Branch("ionNum", &ni);
		tAvalanche->Branch("ionBackNum", &ionBackNum);
		tAvalanche->Branch("ionStartPoints", &ionStartPoints);
		tAvalanche->Branch("ionEndPoints", &ionEndPoints);
		tAvalanche->Branch("ionEndPointsX", &ionEndPointsX);
		tAvalanche->Branch("ionEndPointsY", &ionEndPointsY);
	}

	// Set the signal binning.
	//const int nTimeBins = 10000;
	const double tStep = 10;   //ns
							   //const double timespace = 1./rate*1.e9;    // in ns
	const double timespace = 2.e6;    // in ns // 2ms so that events don't overlap (ion back flow super slow to collect)
									  //const double timespace = 1.e2;
	const double tStart =  0.;
	//const double tEnd = int(nEvents * timespace + ionDelay);
	const double tEnd = int(nEvents * timespace);
	//const double tStep = (tEnd - tStart) / nTimeBins;
	const int nTimeBins = (tEnd - tStart)/tStep;
	sensor->SetTimeWindow(tStart, tStep, nTimeBins);
	//return 0;
	
	// Create an avalanche object
	AvalancheMicroscopic* aval = new AvalancheMicroscopic();
	aval->SetSensor(sensor);
	aval->EnableSignalCalculation();
	
	AvalancheMC* drift = new AvalancheMC();
	if (computeIBF) {
		drift->SetSensor(sensor);
		drift->SetDistanceSteps(2.e-4);
		drift->EnableSignalCalculation();
	}
	
	int division = int(nEvents/10);
	if (nEvents < 10 || useFeSource) division = 1;
	
	//return 0;
	for (int i = 0; i < nEvents; ++i) {
		if (i % division == 0) {
			cout << "\n\n\n\nStart computing event: " << i+1 << "/" << nEvents << endl;
			time_t t = time(NULL);
			PrintTime(t0, t);
			tAvalanche->Write("", TObject::kOverwrite);
		}
		// Initial coordinates of the photon.
		double x0 = width/2. + RndmUniform() * pitch;
		double y0 = depth/2. + RndmUniform() * pitch;
		double z0 = damp + 5*radius + (ddrift-damp-5*radius)*RndmUniform();
		double t0 = i * timespace;
		int ne = 1;	// number of primary electrons created by the photon, if no Fe source then one ionisation only is generated
		if (useFeSource) {
			// Sample the photon energy, using the relative intensities according to XDB.
			const double r = 167. * RndmUniform();
			const double egamma = r < 100. ? 5898.8 : r < 150. ? 5887.6 : 6490.4;
			double dx0 = -1+RndmUniform()*2;	// angle aléatoire du photon
			double dy0 = -1+RndmUniform()*2;
			track.TransportPhoton(x0, y0, z0, t0, egamma, dx0, dy0, -1, ne);
			if (ne < 1) {i--; continue;}	// the detector is quite thin so most of photons go through without converting into electrons, we are not interested in these events
			cout << "number of primary electrons created by the photon = " << ne << endl;
		}
		double e = 0;
		ionBackNum = 0;
		nWinners = 0;
		double x, y, z, t, dx, dy, dz;	// position where the photon creates electrons
		for (int j = 0; j < ne; j++) {	// number of primary electrons created by the photon, = 1 when no photon source
			if (useFeSource) {
				track.GetElectron(j, x, y, z, t, e, dx, dy, dz);
				aval->AvalancheElectron(x, y, z, t, e, 0, 0, -1);
			}
			else aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
			ni = 0;
			neAval = 0;
			aval->GetAvalancheSize(neAval, ni);
			cout << "\nAvalanche size = " << neAval << endl;
			const int np = aval->GetNumberOfElectronEndpoints();
			double xe1, ye1, ze1, te1, e1;
			double xe2, ye2, ze2, te2, e2;
			double xi1, yi1, zi1, ti1;
			double xi2, yi2, zi2, ti2;
			int status;
			for (int j = np; j--;) {	// loop over amplification electrons
				aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
				if (ze2 < 0.008) nWinners++;
				if (computeIBF) {
					electronStartPoints.push_back(ze1);
					electronEndPoints.push_back(ze2);
					drift->DriftIon(xe1, ye1, ze1, te1);
					drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
					if (zi2 > damp+ (ddrift-damp)*0.5) ionBackNum+=1;
					ionStartPoints.push_back(zi1);
					ionEndPoints.push_back(zi2);
					ionEndPointsX.push_back(xi2);
					ionEndPointsY.push_back(yi2);
				}	// end of if (computeIBF)
			} // end of for (int j = np; j--;) (loop over all amplification electrons)
			if (!useFeSource) cout << "nWinners = " << nWinners << " / " << neAval << endl;
			neAvalVec.push_back(neAval);
			nWinnersVec.push_back(nWinners);
		}	// end of for (int j = 0; j < ne; j++) (loop over all primary electrons)
		//cout << ibfRatio << endl;
		tAvalanche->Fill();
		neAvalVec.clear();
		nWinnersVec.clear();
		if (computeIBF) {
			electronStartPoints.clear();
			electronEndPoints.clear();
			ionStartPoints.clear();
			ionEndPoints.clear();
			ionEndPointsX.clear();
			ionEndPointsY.clear();
		}
	}
	tAvalanche->Write("", TObject::kOverwrite);
	
	//cout << "induced charge for electrode " << 2 << " = " << sensor->GetInducedCharge("V2") << endl;
	
	if (computeIBF) {
		for (int k = 0; k < electrodeNum; k++) {
			TTree *tSignal = new TTree(Form("tSignal_%d",k+2),"Currents");
			Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
			Double_t fctInt = 0, fceInt = 0, fciInt = 0;
			if (keepSignal) {
				tSignal->Branch("time", &ft, "time/D");
				tSignal->Branch("totalCurrent", &fct, "totalCurrent/D");
				tSignal->Branch("electronCurrent", &fce, "electronCurrent/D");
				tSignal->Branch("ionCurrent", &fci, "ionCurrent/D");
				tSignal->Branch("totalCurrentInt",&fctInt,"totalCurrentInt/D");     // Integral
				tSignal->Branch("electronCurrentInt",&fceInt,"electronCurrentInt/D");
				tSignal->Branch("ionCurrentInt",&fciInt,"ionCurrentInt/D");
			}
			
			TTree *tCharge = new TTree(Form("tInducedCharge_%d",k+2),"InducedCharge");
			Double_t tic = 0., eic = 0., iic = 0.;
			tCharge->Branch("totalInducedCharge", &tic, "totalInducedCharge/D");
			tCharge->Branch("electronInducedCharge", &eic, "electronInducedCharge/D");
			tCharge->Branch("ionInducedCharge", &iic, "ionInducedCharge/D");
			
			
			//int evNum = 0;
			for (int j = 0; j < nTimeBins; j++) {
				ft = j * tStep;
				//cout << ft << endl;
				fct = sensor->GetSignal(Form("V%d", k+2), j) / ElementaryCharge;
				fce = sensor->GetElectronSignal(Form("V%d", k+2), j) / ElementaryCharge;
				fci = sensor->GetIonSignal(Form("V%d", k+2), j) / ElementaryCharge;
				fctInt+=fct*tStep;
				fceInt+=fce*tStep;
				fciInt+=fci*tStep;
				if (keepSignal) tSignal->Fill();
				
				tic+=fct*tStep;
				eic+=fce*tStep;
				iic+=fci*tStep;
				
				if (int((j+1) * tStep/timespace) == (j+1) * tStep/timespace) {
					//tAvalanche->GetEntry(evNum);
					//cout << "number of amplification electrons = " << nWinners << endl;
					//cout << "total charge induced = " << int(tic) << endl;
					//evNum++;
					tCharge->Fill();
					tic = 0; eic = 0; iic = 0;
				}
			}
			if (keepSignal) tSignal->Write("", TObject::kOverwrite);
			tCharge->Write("", TObject::kOverwrite);
			//cout << "induced charge for electrode " << k+2 << " = " << sensor->GetInducedCharge(Form("V%d", k+2)) << endl;
		}
	}
	f->Close();
	
	
	// Plot the signal.
	const bool plotSignal = false;
	if (plotSignal) {
		TCanvas* cSignal = new TCanvas("signal", "Signal");
		ViewSignal* vSignal = new ViewSignal();
		gPad->SetLeftMargin(0.15);
		vSignal->SetSensor(sensor);
		vSignal->SetCanvas(cSignal);
		//void PlotSignal(const string& label, const bool total = true, const bool electron = false, const bool ion = false);
		vSignal->SetRangeX(0, nEvents*timespace+1000);
		vSignal->PlotSignal("V2", true, true, true);    // signal de la mesh
		cSignal->SaveAs("Figures/Signal.pdf");
	}
	
	
	time_t t1 = time(NULL);
	
	cout << "\n" << nEvents << " events simulated" << endl;
	PrintTime(t0, t1);
	
	//app.Run(true);
}


