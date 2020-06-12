#include <iostream>
#include <fstream>
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

//____________________________________________
void PrintTime(time_t t0, time_t t1)
{
	double tDay = int((t1-t0)/(24*60*60.));
	double tH = int((t1-t0)/(60*60.) - tDay*24);
	double tMin = int((t1-t0)/(60.) - tDay*24*60 - tH*60);
	double tSec = (t1-t0) - tDay*24*60*60 - tH*60*60 - tMin*60;
	std::cout << "Duration of the simulation = " << tDay << " days " << tH << "h " << tMin << " min " << tSec << "s." << std::endl;
	
	std::cout << std::endl;
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
	if (n < 1) {std::cout << "empty folder" << std::endl; return 0;}
	else {
		while (n--) { if (strstr(namelist[n]->d_name, name) != NULL) num++;}
		return num;
	}
}

Int_t GetPrimary(std::string gasName) {
	Int_t nPrimaryTh;
	if (gasName=="Ar-iC4H10") nPrimaryTh = 225;
	else if (gasName=="Ne") nPrimaryTh = 169;       //157 d'apres les calculs...
	else if (gasName=="Ar-CO2") nPrimaryTh = 218;   //222 d'apres les calculs...
	else {std::cout << "What gas??" << std::endl; return 0;}
	return nPrimaryTh;
}

void LoadElectrodeMap(int modelNum, std::map <std::string, int>& electrodeMap) {
	electrodeMap = {};
	if (modelNum == 1 || modelNum == 16 || modelNum == 17 || modelNum == 18) {
		electrodeMap["drift"] = 3;
		electrodeMap["mesh"] = 2;
		electrodeMap["pad"] = 4;
	}
	else if (modelNum >= 2 && modelNum < 5) {   // check it's the same for models 2, 3, 4
		electrodeMap["drift"] = 5;
		electrodeMap["DM up"] = 4;
		electrodeMap["DM down"] = 3;
		electrodeMap["pad"] = 2;
	}
	else if (modelNum >= 5 && modelNum < 8) {   // check it's the same for models 5, 6, 7
		electrodeMap["drift"] = 5;
		electrodeMap["DM up"] = 4;
		electrodeMap["DM down"] = 3;
		electrodeMap["mesh"] = 6;
		electrodeMap["pad"] = 2;
	}
	else if (modelNum >= 8 && modelNum < 10) {   // check it's the same for models 8, 9
		electrodeMap["drift"] = 3;
		electrodeMap["GEM up"] = 6;
		electrodeMap["GEM down"] = 5;
		electrodeMap["mesh"] = 2;
		electrodeMap["pad"] = 4;
	}
	else if (modelNum >= 10 && modelNum < 14) {   // check it's the same for models 10, 11
		electrodeMap["drift"] = 3;
		electrodeMap["mesh top"] = 7;
		electrodeMap["GEM up"] = 6;
		electrodeMap["GEM down"] = 5;
		electrodeMap["mesh"] = 2;
		electrodeMap["pad"] = 4;
	}
	else if (modelNum == 14) {
		electrodeMap["drift"] = 5;
		electrodeMap["DM up"] = 4;
		electrodeMap["DM down"] = 3;
		electrodeMap["pad"] = 2;
	}
	else if (modelNum == 15) {
		electrodeMap["drift"] = 3;
		electrodeMap["GEM2 up"] = 8;
		electrodeMap["GEM2 down"] = 7;
		electrodeMap["GEM1 up"] = 6;
		electrodeMap["GEM1 down"] = 5;
		electrodeMap["mesh"] = 2;
		electrodeMap["pad"] = 4;
	}
	else {std::cout << "no info for this model" << std::endl;}
}

int GetElectrodeNum(int modelNum) {
	std::map <std::string, int> electrodeMap;
	LoadElectrodeMap(modelNum, electrodeMap);
	return electrodeMap.size();
}



bool LoadVariables(int& modelNum, std::string& gasName, int& nEvents, bool& computeIBF)
//T LoadVariable(string elementString)
{
	//std::cout << "elementString = " << elementString << std::endl;
	ifstream file("input.txt", ios::in);
	//T* element = nullptr;
	string a, b;
	bool gasFound = false, modelFound = false, nEventsFound = false, computeInfoFound = false;
	if(file) {
		string line {};
		getline(file, line);	//first line does not contains info
		while(getline(file, line)) {
			if (line.find("modelNum") != std::string::npos) {
				stringstream stream(line);
				stream >> a >> b >> modelNum;
				modelFound = true;
			}
			else if (line.find("gasName") != std::string::npos) {
				stringstream stream(line);
				stream >> a >> b >> gasName;
				gasFound = true;
			}
			else if (line.find("nEvents") != std::string::npos) {
				stringstream stream(line);
				stream >> a >> b >> nEvents;
				nEventsFound = true;
			}
			else if (line.find("computeIBF") != std::string::npos) {
				stringstream stream(line);
				stream >> a >> b >> computeIBF;
				computeInfoFound = true;
			}
		}
		//std::cout << "Element " << elementString << " not found in input.txt";
		if (modelFound && gasFound && nEventsFound && computeInfoFound) {
			std::cout << std::endl;
			std::cout << "#################################" << std::endl;
			std::cout << "#\tmodelNum = " << modelNum << "\t\t#" << std::endl;
			if (gasName.length() < 5) std::cout << "#\tgasName = " << gasName << "\t\t#" << std::endl;
			else std::cout << "#\tgasName = " << gasName << "\t#" << std::endl;
			std::cout << "#\tnEvents = " << nEvents << "\t\t#" << std::endl;
			std::cout << "#\tcomputeIBF = " << computeIBF << "\t\t#" << std::endl;
			std::cout << "#################################" << std::endl;
			std::cout << std::endl;
			return true;
		}
		else return false;
	}
	else cout << "Error: not possible to open input.txt file in reading mode" << endl;
	return false;
	//return *element;
}

bool LoadVariables(int& modelNum, std::string& gasName) {
	int nEvents = 0;
	bool computeIBF = false;
	return LoadVariables(modelNum, gasName, nEvents, computeIBF);
}

int GetMaxModelNum() {return 18;}
