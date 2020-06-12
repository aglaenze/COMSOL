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

// Set up detector geometry

void LoadParameters(int modelNum, double& damp, double& ddrift, double& radius, double& pitch, double& width, double& depth, std::vector<double>& zElectrodes) {
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
		zElectrodes = {ddrift, damp, 0.0128, 0};
    }
    else if (modelNum == 9) {   // MM + GEM
        damp = 0.0128+0.4+0.0060;
        radius = 0.0009;   // cm
        ddrift = 0.7;  //cm
        pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128, 0};
    }
    else if (modelNum == 10) {  // MM + (GEM+MM) combined
        damp = 0.0128+0.2+0.0060+0.0128;
        radius = 0.0009;   // cm
        ddrift = 0.5;  //cm
        pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.2+0.0060, 0.0128, 0};
    }
    else if (modelNum == 11) {   // MM + (GEM+MM) combined
        damp = 0.0128+0.4+0.0060+0.0128;
        radius = 0.0009;   // cm
        ddrift = 0.7;  //cm
        pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.4+0.0060, 0.0128, 0};
    }
	else if (modelNum == 12) {   // MM + (GEM+MM) combined
		damp = 0.0128+0.2+0.0060+0.0256;
		radius = 0.0009;   // cm
		ddrift = 0.7;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.2+0.0060, 0.0128, 0};
	}
	else if (modelNum == 13) {   // MM + (GEM+MM) combined
		damp = 0.0128+0.4+0.0060+0.0256;
		radius = 0.0009;   // cm
		ddrift = 0.7;  //cm
		pitch = 0.0126;    // cm
		zElectrodes = {ddrift, damp, 0.0128+0.4+0.0060, 0.0128, 0};
	}
	else if (modelNum == 14) {   // DM
		damp = 0.0540;
		radius = 0.0009;   // cm
		ddrift = 0.35;  //cm
		pitch = 0.0040;    // cm
		zElectrodes = {ddrift, damp, 0.0220, 0};
	}
	else if (modelNum == 15) {   // MM + 2 GEMs
		damp = 0.6125;	// cm
		radius = 0.0009;   // cm
		ddrift = 1.4125;  //cm
		pitch = 0.0180;    // cm
		zElectrodes = {ddrift, damp, 0.0125+0.4, 0.0125, 0};
	}
    else {std::cout << "Model num?" << std::endl; return;}
    width = periodicityNum * pitch;
    depth = periodicityNum * pitch;
}

void LoadParameters(int modelNum, double& damp, double& ddrift, double& radius, double& pitch, double& width, double& depth) {
	std::vector<double> zElectrodes = {};
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth, zElectrodes);
}

void LoadParameters(int modelNum, std::vector<double>& zElectrodes) {
	double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
	LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth, zElectrodes);
}

void DrawDetector(int modelNum, std::vector<int> hvList) {
	// Draw geometry with voltages
	
	gStyle->SetTextSize(0.06);
	
	std::map <std::string, int> electrodeMap;
	LoadElectrodeMap(modelNum, electrodeMap);
	int electrodeNum = GetElectrodeNum(modelNum);
	
	// insert pad voltage in HV list
	hvList.insert(hvList.begin(),0);
	
	double yMax = 0.8; // coord of drift electrode
	double space = yMax/(electrodeNum-1);
	
	TBox* detectorBox = new TBox(0,0,1,1);
	detectorBox->Draw();
	
	int i = 0;	// i is the stage index
	
	// draw the electrodes from top to bottom (as declared in the function LoadElectrodeMap)
	
	std::map<std::string, int>::iterator it = electrodeMap.begin();
	// Iterate over the map using Iterator until end.
	while (it != electrodeMap.end()) {
		double yCoord = yMax-i*space;
		TLine* electrodeLine = new TLine(0, yCoord, 1, yCoord);
		electrodeLine->SetLineStyle(9);
		if (i == 0 || i == electrodeNum-1) electrodeLine->SetLineStyle(1);
		TLatex* electrodeText = new TLatex(0.1, yCoord+0.03, Form("V_{%s} = %d V", (it->first).c_str(), hvList[electrodeNum-1-i]));
		electrodeLine->Draw("same");
		electrodeText->Draw("same");
		// Increment the Iterator to point to next entry
		it++;
		i++;
	}
	
}

