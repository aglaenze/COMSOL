#include <iostream>
#include <fstream>
#include <cmath>

#include <string.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TString.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

using namespace std;

/* Set up detector geometry */

void LoadParameters(int modelNum, double& damp, double& ddrift, double& radius, double& pitch, double& width, double& depth, vector<double>& zElectrodes) {
	int periodicityNum = 2000;
	if (modelNum == 1 || modelNum == 16 || modelNum == 17 || modelNum == 18) {
		damp = 0.0128;
		ddrift = 0.5;      // cm
		radius = 0.0009;   // cm
		pitch = 0.0063;    // cm
		if (modelNum == 16) pitch = 0.0078;
		else if (modelNum == 17) pitch = 0.0096;
		else if (modelNum == 18) pitch = 0.0045;
		zElectrodes = {ddrift, damp, 0};
	}
	else if (modelNum == 2 || modelNum == 3 || modelNum == 4) {
		damp = 2*0.0128;
		radius = 0.0009;   // cm
		ddrift = 0.5;  //cm
		if (modelNum == 2) pitch = 0.0045;    // cm
		else if (modelNum == 3) pitch = 0.0078;    // cm
		else if (modelNum == 4) pitch = 0.0063;    // cm
		zElectrodes = {ddrift, damp, 0.0128, 0};
	}
	else if (modelNum >= 5 && modelNum < 8) {
		damp = 2*0.0128+0.2;
		radius = 0.0009;   // cm
		ddrift = 0.7;  //cm
		pitch = 0.0063;    // cm
		zElectrodes = {ddrift, damp, 0.2+0.0128, 0.0128, 0};
	}
	else if (modelNum == 8) {  // MM + GEM
		damp = 0.0128+0.2+0.0060;
		radius = 0.0009;   // cm
		ddrift = 0.5;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.2, 0.0128, 0};
	}
	else if (modelNum == 9) {   // MM + GEM
		damp = 0.0128+0.4+0.0060;
		radius = 0.0009;   // cm
		ddrift = 0.7;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.4, 0.0128, 0};
	}
	else if (modelNum == 10 || modelNum == 19) {  // MM + (GEM+MM) combined
		damp = 0.0128+0.2+0.0060+0.0128;
		if (modelNum == 19) damp = 0.0128+0.2+0.0060+0.0064;
		radius = 0.0009;   // cm
		ddrift = 0.5;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.2+0.0060, 0.0128+0.2, 0.0128, 0};
	}
	else if (modelNum == 11) {   // MM + (GEM+MM) combined
		damp = 0.0128+0.4+0.0060+0.0128;
		radius = 0.0009;   // cm
		ddrift = 0.7;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.4+0.0060, 0.0128+0.4, 0.0128, 0};
	}
	else if (modelNum == 12) {   // MM + (GEM+MM) combined
		damp = 0.0128+0.2+0.0060+0.0256;
		radius = 0.0009;   // cm
		ddrift = 0.5;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.2+0.0060, 0.0128+0.2, 0.0128, 0};
	}
	else if (modelNum == 13) {   // MM + (GEM+MM) combined
		damp = 0.0128+0.4+0.0060+0.0256;
		radius = 0.0009;   // cm
		ddrift = 0.7;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.4+0.0060, 0.0128+0.4, 0.0128, 0};
	}
	else if (modelNum == 14) {   // DM
		damp = 0.0540;
		radius = 0.0009;   // cm
		ddrift = 0.35;  //cm
		pitch = 0.0040;    // cm
		zElectrodes = {ddrift, damp, 0.0220, 0};
	}
	else if (modelNum == 15) {   // MM + 2 GEMs
		damp = 0.6125+2*0.0070;	// cm
		radius = 0.0009;   // cm
		ddrift = 1.4125;  //cm
		pitch = 0.0180;    // cm
		zElectrodes = {ddrift, damp, 0.6125+0.0070, 0.4125+0.0070, 0.0125+0.4, 0.0125, 0};
	}
	else if (modelNum == 20) {   // MGEM1
		damp = 0.0128+0.4+0.0060+0.0660;
		radius = 0.0009;   // cm
		ddrift = 1.0;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.4+0.0060, 0.0128+0.4, 0.0128, 0};
	}
	else if (modelNum == 21) {   // MGEM3
		damp = 0.0128+0.4+0.0060+0.0128;
		radius = 0.0013;   // cm
		ddrift = 0.7;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.4+0.0060, 0.0128+0.4, 0.0128, 0};
	}
	else {std::cout << "Model num?" << std::endl; return;}
	width = periodicityNum * pitch;
	depth = periodicityNum * pitch;
}

void LoadParameters(int modelNum, double& damp, double& ddrift, double& radius, double& pitch, double& width, double& depth) {
	vector<double> zElectrodes = {};
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth, zElectrodes);
}

void LoadParameters(int modelNum, vector<double>& zElectrodes) {
	double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth, zElectrodes);
}

void DrawDetector(int modelNum, vector<int> hvList) {
	// Draw geometry with voltages
	
	gStyle->SetTextSize(0.05);
	
	map <string, int, NoSorting> electrodeMap;
	LoadElectrodeMap(modelNum, electrodeMap);
	int electrodeNum = GetElectrodeNum(modelNum);
	
	vector<double> zElectrodes = {};
	LoadParameters(modelNum, zElectrodes);
	
	// insert pad voltage in HV list
	hvList.insert(hvList.begin(),0);
	
	double yMax = 0.8; // coord of drift electrode
	double space = yMax/(electrodeNum-1);
	
	TBox* detectorBox = new TBox(0,0,1,1);
	detectorBox->Draw();
	
	int i = 0;	// i is the stage index
	
	// draw the electrodes from top to bottom (as declared in the function LoadElectrodeMap)
	
	map<string, int>::iterator it = electrodeMap.begin();
	// Iterate over the map using Iterator until end.
	while (it != electrodeMap.end()) {
		double yCoord = yMax-i*space;
		TLine* electrodeLine = new TLine(0, yCoord, 1, yCoord);
		electrodeLine->SetLineStyle(9);
		// 0 = drift
		if (i == 0 || i == electrodeNum-1) electrodeLine->SetLineStyle(1);
		TLatex* electrodeText = new TLatex(0.05, yCoord+0.03, Form("V_{%s} = %d V", (it->first).c_str(), hvList[electrodeNum-1-i]));
		TText* zElectrode = new TText(0.6, yCoord+0.03, Form("z = %.3f mm", zElectrodes[i]*10));
		electrodeLine->Draw("same");
		electrodeText->Draw("same");
		zElectrode->Draw("same");
		// Computation of electric field
		if (i != electrodeNum-1) {
			double electricField = (hvList[electrodeNum-1-i]-hvList[electrodeNum-1-(i+1)])/(zElectrodes[i]-zElectrodes[i+1]);
			TString fieldtxt = Form("E = %.1f V/cm", electricField);
			if (electricField>1000.) fieldtxt = Form("E = %.2f kV/cm", electricField/1000.);
			TText* electricFieldtxt = new TText(0.05, yCoord-0.5*space, fieldtxt);
			electricFieldtxt->SetTextColor(kBlue);
			electricFieldtxt->Draw("same");
		}
		// Increment the Iterator to point to next entry
		it++;
		i++;
	}
	
}

void DrawElectrodes(int modelNum, double zMin, double zMax, bool yDir = false) {
	map <string, int, NoSorting> electrodeMap;
	LoadElectrodeMap(modelNum, electrodeMap);
	vector<double> zElectrodes = {};
	LoadParameters(modelNum, zElectrodes);
	
	double zPadMin = gPad->GetUymin();
	double zPadMax = gPad->GetUymax();

	double scale = (zPadMax-zPadMin)/(zMax-zMin);

	for (int i = 0; i<(int)zElectrodes.size(); i++) {
		TLine* electrodeLine;
		double newZ = zPadMin+ (zElectrodes[i]-zMin)*scale;
		if (yDir) {electrodeLine = new TLine(zElectrodes[i], zMin, zElectrodes[i], zMax);}
		else {electrodeLine = new TLine(gPad->GetUxmin(), newZ, gPad->GetUxmax(), newZ);}
		//std::cout << zElectrodes[i] << std::endl;
		electrodeLine->SetLineStyle(9);
		if (yDir) {electrodeLine->SetLineStyle(3); electrodeLine->SetLineColor(15);}
		electrodeLine->Draw("same");
	}
	
}
