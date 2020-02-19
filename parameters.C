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

void LoadParameters(int modelNum, int& periodicityNum, double& damp, double& ddrift, double& dmylar, double& radius, double& pitch, double& width, double& depth) {
    const double dPA = 0.024 + 0.0027;
    const double dSA = 0.012;
    periodicityNum = 2000;
    if (modelNum == 1) {
        damp = 0.0128;
        ddrift = 0.5;      // cm
        radius = 0.0009;   // cm
        pitch = 0.0063;    // cm
    }
    else if (modelNum == 4) {
        damp = dSA + dPA;
        radius = 0.0027;   // cm
        ddrift = 0.3 + damp + radius;  //cm
        pitch = 0.025;    // cm
    }
    else if (modelNum == 5) {
        damp = 0.0120;
        radius = 0.0009;   // cm
        ddrift = 0.5182;  //cm
        pitch = 0.0063;    // cm
    }
    else if (modelNum == 6) {
        damp = 0.0120;
        radius = 0.0009;   // cm
        ddrift = 0.5240;  //cm
        pitch = 0.0063;    // cm
    }
    else if (modelNum == 7) {
        damp = 2*0.0128;
        radius = 0.0009;   // cm
        ddrift = 0.5;  //cm
        pitch = 0.0045;    // cm
    }
    else if (modelNum == 8) {
        damp = 2*0.0128;
        radius = 0.0009;   // cm
        ddrift = 0.5;  //cm
        pitch = 0.0078;    // cm
    }
    else if (modelNum == 9) {
        damp = 2*0.0128;
        radius = 0.0009;   // cm
        ddrift = 0.5;  //cm
        pitch = 0.0063;    // cm
    }
    else if (modelNum >= 10 && modelNum < 13) {
        damp = 2*0.0128+0.2;
        radius = 0.0009;   // cm
        ddrift = 0.7;  //cm
        pitch = 0.0063;    // cm
    }
    else {std::cout << "Model num?" << std::endl; return;}
    dmylar = 3.;       // cm
    width = periodicityNum * pitch;
    depth = periodicityNum * pitch;
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
        //gas->SetComposition("Ar", 93., "CO2", 7.);
        gas->SetComposition("Ar", 90., "CO2", 10.);
        //rPenning = 0.5;
        rPenning = 0.4;
        gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
        gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
        //gas->LoadIonMobility(path + "/Data/IonMobility_CO2+_CO2.txt");
    }
    else {std::cout << "What gas??" << std::endl; return 0;}
    gas->SetTemperature(293.15);
    gas->SetPressure(Garfield::AtmosphericPressure);
    gas->EnableDrift();
    //gas->EnablePenningTransfer(rPenning, lambdaPenning);
    return gas;
}

Garfield::ComponentComsol* InitiateField(int modelNum, int hvMesh, int hvDrift, Garfield::MediumMagboltz* gas) {
    // Load the field map.
    std::string dataFolder = Form("COMSOL_data/model%d/", modelNum);
    std::string dataFile = dataFolder + Form("ewfield-%d-%d.txt", hvMesh, hvDrift);
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

Garfield::ComponentComsol* InitiateField(int modelNum, int hvDmDown, int hvDmUp, int hvDrift, Garfield::MediumMagboltz* gas) {
    // Load the field map.
    std::string dataFolder = Form("COMSOL_data/model%d/", modelNum);
    std::string dataFile = dataFolder + Form("ewfield-%d-%d-%d.txt", hvDmDown, hvDmUp, hvDrift);
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

Garfield::ComponentComsol* InitiateField(int modelNum, int hvMesh, int hvDmDown, int hvDmUp, int hvDrift, Garfield::MediumMagboltz* gas) {
    // Load the field map.
    std::string dataFolder = Form("COMSOL_data/model%d/", modelNum);
    std::string dataFile = dataFolder + Form("ewfield-%d-%d-%d-%d.txt", hvMesh, hvDmDown, hvDmUp, hvDrift);
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

