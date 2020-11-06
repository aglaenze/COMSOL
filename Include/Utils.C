#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <dirent.h>
#include <string.h>
#include <stdio.h>

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

struct NoSorting {
	bool operator()(string const& lhsIn, string const& rhsIn) const {
		if (lhsIn > rhsIn) return true;	// if condition to silence the warning
		else return true;
	}
};

//____________________________________________
void PrintTime(time_t t0, time_t t1)
{
	double tDay = int((t1-t0)/(24*60*60.));
	double tH = int((t1-t0)/(60*60.) - tDay*24);
	double tMin = int((t1-t0)/(60.) - tDay*24*60 - tH*60);
	double tSec = (t1-t0) - tDay*24*60*60 - tH*60*60 - tMin*60;
	cout << "Duration of the simulation = " << tDay << " days " << tH << "h " << tMin << " min " << tSec << "s." << endl;
	
	cout << endl;
}

TLatex* PutText( Double_t x_ndc, Double_t y_ndc, TString value, Double_t fontSize = 0.04 )
{
	TLatex* text = new TLatex();
	text->SetNDC( true );
	text->SetTextSize( fontSize );
	text->DrawLatex( x_ndc, y_ndc, value );
	return text;
}

int GetNumberOfFiles(TString path, TString name) {
	Int_t num = 0;
	struct dirent **namelist;
	Int_t n = scandir(path, &namelist, 0, alphasort);
	if (n < 1) {cout << "empty folder" << endl; return 0;}
	else {
		while (n--) { if (strstr(namelist[n]->d_name, name) != NULL) num++;}
		return num;
	}
}

Int_t GetPrimary(string gasName) {
	Int_t nPrimaryTh;
	if (gasName=="Ar-iC4H10") nPrimaryTh = 225;
	else if (gasName=="Ne") nPrimaryTh = 169;       //157 d'apres les calculs...
	else if (gasName=="Ar-CO2") nPrimaryTh = 218;   //222 d'apres les calculs...
	else {cout << "What gas??" << endl; return 0;}
	return nPrimaryTh;
}

void LoadElectrodeMap(int modelNum, map <string, int, NoSorting>& electrodeMap) {
	electrodeMap = {};
	if (modelNum == 1 || modelNum == 16 || modelNum == 17 || modelNum == 18) {
		electrodeMap["pad"] = 4;
		electrodeMap["mesh"] = 2;
		electrodeMap["drift"] = 3;
	}
	else if (modelNum >= 2 && modelNum < 5) {
		electrodeMap["pad"] = 2;
		electrodeMap["DM down"] = 3;
		electrodeMap["DM up"] = 4;
		electrodeMap["drift"] = 5;
	}
	else if (modelNum >= 5 && modelNum < 8) {
		electrodeMap["pad"] = 2;
		electrodeMap["mesh"] = 6;
		electrodeMap["DM down"] = 3;
		electrodeMap["DM up"] = 4;
		electrodeMap["drift"] = 5;
	}
	else if (modelNum >= 8 && modelNum < 10) {
		electrodeMap["pad"] = 4;
		electrodeMap["mesh"] = 2;
		electrodeMap["GEM down"] = 5;
		electrodeMap["GEM up"] = 6;
		electrodeMap["drift"] = 3;
	}
	else if ((modelNum >= 10 && modelNum < 14) || (modelNum >= 19 && modelNum < 22)) {
		electrodeMap["pad"] = 4;
		electrodeMap["mesh"] = 2;
		electrodeMap["GEM down"] = 5;
		electrodeMap["GEM up"] = 6;
		electrodeMap["mesh top"] = 7;
		electrodeMap["drift"] = 3;
	}
	else if (modelNum == 14) {
		electrodeMap["pad"] = 2;
		electrodeMap["DM down"] = 3;
		electrodeMap["DM up"] = 4;
		electrodeMap["drift"] = 5;
	}
	else if (modelNum == 15) {
		electrodeMap["pad"] = 4;
		electrodeMap["mesh"] = 2;
		electrodeMap["GEM1 down"] = 5;
		electrodeMap["GEM1 up"] = 6;
		electrodeMap["GEM2 down"] = 7;
		electrodeMap["GEM2 up"] = 8;
		electrodeMap["drift"] = 3;
	}
	else {cout << "No info for this model in LoadElectrodeMap" << endl;}
}

int GetElectrodeNum(int modelNum) {
	map <string, int, NoSorting> electrodeMap;
	LoadElectrodeMap(modelNum, electrodeMap);
	return electrodeMap.size();
}

bool LoadVariables(int& modelNum, string& gasName, int& nEvents, bool& computeIBF, bool& useFeSource, bool& plotDrift2D, bool& plotDrift3D)
//T LoadVariable(string elementString)
{
	//cout << "elementString = " << elementString << endl;
	ifstream file("input.txt", ios::in);
	//T* element = nullptr;
	string a, b;
	bool gasFound = false, modelFound = false, nEventsFound = false, computeInfoFound = false, feInfoFound = false, plot2dInfoFound = false, plot3dInfoFound = false;
	if(file) {
		string line {};
		getline(file, line);	//first line does not contains info
		while(getline(file, line)) {
			if (line.find("modelNum") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> modelNum;
				modelFound = true;
			}
			else if (line.find("gasName") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> gasName;
				gasFound = true;
			}
			else if (line.find("nEvents") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> nEvents;
				nEventsFound = true;
			}
			else if (line.find("computeIBF") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> computeIBF;
				computeInfoFound = true;
			}
			else if (line.find("useFeSource") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> useFeSource;
				feInfoFound = true;
			}
			else if (line.find("plotDrift2D") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> plotDrift2D;
				plot2dInfoFound = true;
			}
			else if (line.find("plotDrift3D") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> plotDrift3D;
				plot3dInfoFound = true;
			}
		}
		//cout << "Element " << elementString << " not found in input.txt";
		return (modelFound && gasFound && feInfoFound && nEventsFound && computeInfoFound && plot2dInfoFound && plot3dInfoFound);
	}
	else cout << "Error: not possible to open input.txt file in reading mode" << endl;
	return false;
	//return *element;
}

bool LoadVariables(int& modelNum, string& gasName) {
	int nEvents = 0;
	bool computeIBF = false, useFeSource = false;
	bool plotDrift2D = 0, plotDrift3D = 0;
	bool result = LoadVariables(modelNum, gasName, nEvents, computeIBF, useFeSource, plotDrift2D, plotDrift3D);
	if (!result) return false;
	cout << endl;
	cout << "#################################" << endl;
	cout << "#\tmodelNum = " << modelNum << "\t\t#" << endl;
	if (gasName.length() < 5) cout << "#\tgasName = " << gasName << "\t\t#" << endl;
	else cout << "#\tgasName = " << gasName << "\t#" << endl;
	cout << "#################################" << endl;
	cout << endl;
	return true;
}

bool LoadVariables(int& modelNum, string& gasName, bool& plotDrift2D, bool& plotDrift3D) {
	int nEvents = 0;
	bool computeIBF = false, useFeSource = false;
	bool result =  LoadVariables(modelNum, gasName, nEvents, computeIBF, useFeSource, plotDrift2D, plotDrift3D);
	if (!result) return false;
	cout << endl;
	cout << "#################################" << endl;
	cout << "#\tmodelNum = " << modelNum << "\t\t#" << endl;
	if (gasName.length() < 5) cout << "#\tgasName = " << gasName << "\t\t#" << endl;
	else cout << "#\tgasName = " << gasName << "\t#" << endl;
	cout << "#\tplotDrift2D = " << plotDrift2D << "\t\t#" << endl;
	cout << "#\tplotDrift3D = " << plotDrift3D << "\t\t#" << endl;
	cout << "#################################" << endl;
	cout << endl;
	return true;
}

bool LoadVariables(int& modelNum, string& gasName, int& nEvents, bool& computeIBF, bool& useFeSource) {
	bool plotDrift2D = 0, plotDrift3D = 0;
	bool result =  LoadVariables(modelNum, gasName, nEvents, computeIBF, useFeSource, plotDrift2D, plotDrift3D);
	if (!result) return false;
	cout << endl;
	cout << "#################################" << endl;
	cout << "#\tmodelNum = " << modelNum << "\t\t#" << endl;
	if (gasName.length() < 5) cout << "#\tgasName = " << gasName << "\t\t#" << endl;
	else cout << "#\tgasName = " << gasName << "\t#" << endl;
	cout << "#\tnEvents = " << nEvents << "\t\t#" << endl;
	cout << "#\tcomputeIBF = " << computeIBF << "\t\t#" << endl;
	cout << "#\tuseFeSource = " << useFeSource << "\t\t#" << endl;
	cout << "#################################" << endl;
	cout << endl;
	return true;
}

int GetMaxModelNum() {
	int num = 0;
	bool loop = true;
	while (loop) {
		num ++;
		string folder = "COMSOL_data/model" + to_string(num);
		const char* path = &folder[0];
		DIR* rep = NULL;
		rep = opendir(path);
		if (rep == NULL) {loop = false; break;}
	}
	return num-1;
}