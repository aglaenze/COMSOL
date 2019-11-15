#include <iostream>
#include <fstream>
#include <cmath>
#include <dirent.h>

#include <string.h>
#include <TLatex.h>
#include <TString.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"


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

