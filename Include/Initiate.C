#include <iostream>
#include <fstream>
#include <cmath>

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/MediumGas.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/ComponentComsol.hh"

using namespace std;
using namespace Garfield;

// Make a gas medium.
MediumMagboltz* InitiateGas(string gasName) {
	
	MediumMagboltz* gas = new MediumMagboltz();
	double rPenning;
	const double lambdaPenning = 0.;    // parameter for sampling the distance of the Penning electron with respect to the excitation
	const string path = getenv("GARFIELD_HOME");
	if (gasName=="Ar-iC4H10") {
		gas->SetComposition("Ar", 95., "C4H10", 5.);
        //rPenning = 0.4753;
        //rPenning = 0.478;
        //rPenning = 0.48;
        rPenning = 0.18;
		gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
		gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
	}
	else if (gasName=="Ne") {
		gas->SetComposition("Ne", 90., "CF4", 10.);
		rPenning = 0.6;
		gas->EnablePenningTransfer(rPenning, lambdaPenning, "ne");
		gas->LoadIonMobility(path + "/Data/IonMobility_Ne+_Ne.txt");
	}
	else if (gasName=="Ar-CO2") {
		gas->SetComposition("Ar", 93., "CO2", 7.);
		//gas->SetComposition("Ar", 90., "CO2", 10.);
		rPenning = 0.6;
		//rPenning = 0.4;
		gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
		gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
		//gas->LoadIonMobility(path + "/Data/IonMobility_CO2+_CO2.txt");
	}
	else {cout << "What gas??" << endl; return 0;}
	gas->SetTemperature(293.15);
	gas->SetPressure(AtmosphericPressure);
	gas->EnableDrift();
	//gas->EnablePenningTransfer(rPenning, lambdaPenning);
	return gas;
}


ComponentComsol* InitiateField(int modelNum, vector<int> hvList, MediumMagboltz* gas, bool remote = false) {
	// Load the field map.
	string dataFolder = Form("COMSOL_data/model%d/", modelNum);
	if (remote) dataFolder = "./";
	string dataFile = dataFolder + "ewfield";
	for (int k = 0; k < (int)hvList.size(); k++) dataFile += Form("-%d", hvList[k]);
	dataFile += ".txt";
	ifstream meshFile(dataFolder+"mesh.mphtxt");
	ifstream materialsFile(dataFolder+"dielectrics.dat");
	ifstream potentialFile(dataFile);
	if (!meshFile) cout << "no mesh file to read" << endl;
	if (!materialsFile) cout << "no materials file to read" << endl;
	if (!potentialFile) cout << "no data file (with V and E) to read" << endl;
	
	if (!(meshFile && materialsFile && potentialFile)) {return nullptr;}
	else {
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
		return fm;
	}
}
