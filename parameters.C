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
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/ComponentComsol.hh"


// Set up detector geometry
const double pitch = 0.0025;    // cm
const double damp = 0.0128;
const double ddrift = 0.5;      // cm
const double dmylar = 3.;       // cm
const double radius = 0.0004;   // cm
const int periodicityNum = 5000;
double width = periodicityNum * pitch;
double depth = periodicityNum * pitch;


Int_t GetPrimary(std::string gasName) {
    Int_t nPrimaryTh;
    if (gasName=="Ar-iC4H10") nPrimaryTh = 225;
    else if (gasName=="Ne") nPrimaryTh = 169;       //157 d'apres les calculs...
    else if (gasName=="Ar-CO2") nPrimaryTh = 218;   //222 d'apres les calculs...
    else {std::cout << "What gas??" << std::endl; return 0;}
    return nPrimaryTh;
}

 // Make a gas medium.
Garfield::MediumMagboltz* InitiateGas(std::string gasName) {
    
    Garfield::MediumMagboltz* gas = new Garfield::MediumMagboltz();
    double rPenning;
    const double lambdaPenning = 0.;    // parameter for sampling the distance of the Penning electron with respect to the excitation
    const std::string path = getenv("GARFIELD_HOME");
    if (gasName=="Ar-iC4H10") {
        gas->SetComposition("Ar", 95., "iC4H10", 5.);
        rPenning = 0.45;
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
        rPenning = 0.4;
        gas->EnablePenningTransfer(rPenning, lambdaPenning, "ne");
        gas->LoadIonMobility(path + "/Data/IonMobility_Ne+_Ne.txt");
    }
    else {std::cout << "What gas??" << std::endl; return 0;}
    gas->SetTemperature(293.15);
    gas->SetPressure(Garfield::AtmosphericPressure);
    gas->EnableDrift();
    //gas->EnablePenningTransfer(rPenning, lambdaPenning);
    return gas;
}

Garfield::ComponentComsol* InitiateField(int modelNum, int hvMesh, Garfield::MediumMagboltz* gas) {
    // Load the field map.
    std::string dataFolder = Form("COMSOL_data/model%d/", modelNum);
    std::string dataFile = dataFolder + Form("ewfield_%dV.txt", hvMesh);
    // Load the field map.
    Garfield::ComponentComsol* fm = new Garfield::ComponentComsol();
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


