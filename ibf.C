#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>

#include "_Utils.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {
    
    //______________________
    // variables
    std::string gasName = "Ar"; // Ar or Ne
    const int modelNum = 1;
    //____________________
    
    time_t t0 = time(NULL);
    
    if (argc < 2) {
        std::cout << "Please enter HVmesh like this: ./ibf $hvMesh" << std::endl;
        return 0;
    }
    const int hvMesh = atoi(argv[1]);
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    // Set up detector geometry
    //GeometrySimple* geo = new GeometrySimple();
    const double pitch = 0.0025;    // cm
    const double damp = 0.0128;
    const double ddrift = 0.5;      // cm
    const double radius = 0.0004;   // cm
    const int periodicityNum = 5000;
    double width = periodicityNum * pitch;
    double depth = periodicityNum * pitch;
    
    
    // Make a gas medium.
    MediumMagboltz* gas = new MediumMagboltz();
    double rPenning;
    const double lambdaPenning = 0.;    // parameter for sampling the distance of the Penning electron with respect to the excitation
    const std::string path = getenv("GARFIELD_HOME");
    if (gasName=="Ar") {
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
    else {std::cout << "What gas??" << std::endl; return 0;}
    //return 0;
    gas->SetTemperature(293.15);
    gas->SetPressure(AtmosphericPressure);
    gas->EnableDrift();

    //gas->EnablePenningTransfer(rPenning, lambdaPenning);
    
    Medium* wire = new Medium();
    wire->SetTemperature(293.15);
    wire->SetPressure(AtmosphericPressure);
    //wire->SetComposition("Stainless steel");
    
    
    // Load the field map.
    std::string dataFolder = Form("COMSOL_data/model%d/", modelNum);
    std::string dataFile = dataFolder + Form("ewfield_%dV.txt", hvMesh);
    // Load the field map.
    ComponentComsol* fm = new ComponentComsol();
    fm->Initialise(dataFolder+"mesh.mphtxt", dataFolder+"dielectrics.dat", dataFile);
    fm->PrintMaterials();
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
    fm->PrintRange();
    
    // Associate the gas with the corresponding field map material.
    
     const unsigned int nMaterials = fm->GetNumberOfMaterials();
     for (unsigned int i = 0; i < nMaterials; ++i) {
     const double eps = fm->GetPermittivity(i);
     if (eps > 1.) fm->SetMedium(i, wire);
     }
    fm->SetMedium(0, gas);
    
    fm->PrintMaterials();
    
    // Make a sensor.
    Sensor sensor;
    sensor.AddComponent(fm);
    sensor.SetArea(0, 0, 0, width, depth, ddrift);
    //sensor.SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
    
    /*
    // To look at the avalanche
    ViewDrift* driftView = new ViewDrift();
    driftView->SetArea(390*pitch, 390*pitch, 0, 410*pitch, 410*pitch, damp+6*pitch);
    */
    
    // Create an avalanche object
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(&sensor);
    AvalancheMC* drift = new AvalancheMC();
    drift->SetSensor(&sensor);
    drift->SetDistanceSteps(2.e-4);
    //aval->EnablePlotting(driftView);
    //drift->EnablePlotting(driftView);
    
    double totalIonNum = 0;
    double totalIonBackNum = 0;
    
    // Create ROOT histograms of the signal and a file in which to store them.
    TH1F* hIbf = new TH1F("hIbf", "hIbf", 100, 0, 10);
    //TH1::StatOverflows(true);
    hIbf->SetXTitle("IBF (%)");
        hIbf->SetYTitle("# counts");
    
    //TH1F* hze1 = new TH1F("hze1", "hze1", 10000, 0, damp);

    // Write the histograms to the TFile.
    const char* name = Form("rootFiles/%s/model%d/ibf_%dV.root", gasName.c_str(), modelNum, hvMesh);
    //const char* name = Form("rootFiles/%s/model%d/test.root", gasName.c_str(), modelNum);
    TFile* f = new TFile(name, "RECREATE");
    
    const int nEvents = 20;
    int division = int(nEvents/20);
    
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (i % 100 == 0) std::cout << "\n\n\n\n" << i << "/" << nEvents << "\n";
        // Initial coordinates of the photon.
        double x0 = width/2. + RndmUniform() * pitch;
        //double y0 = RndmUniform() * depth;
        double y0 = depth/2. + RndmUniform() * pitch;
        double z0 = damp + radius + (ddrift-damp-radius)*RndmUniform();
        double t0 = 0;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        int ne2 = 0, ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        std::cout << "Avalanche size = " << ne2 << std::endl;
        if (ne2 < 2) continue;
        //hElectrons->Fill(ne2);
        // to look at avalanche
            const int np = aval->GetNumberOfElectronEndpoints();
            double xe1, ye1, ze1, te1, e1;
            double xe2, ye2, ze2, te2, e2;
            double xi1, yi1, zi1, ti1;
            double xi2, yi2, zi2, ti2;
            int status;
            double ionBackNum = 0;
            for (int j = np; j--;) {
                aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
                //std::cout << "departure of the electron in x y z : " << xe1 << " " << ye1 << " " <<  ze1 << std::endl;
                //std::cout << "arrival of the electron in x y z : " << xe2 << " " << ye2 << " " <<  ze2 << std::endl;
                drift->DriftIon(xe1, ye1, ze1, te1);
                //hze1->Fill(ze1);
                drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
                //std::cout << "arrival of the ion in z : " << zi2 << std::endl;
                if (zi2 > 3*damp) ionBackNum+=1;
            }
            hIbf->Fill( ionBackNum*100./ni);
        if (i % division == 0) hIbf->Write("", TObject::kOverwrite);
            totalIonBackNum += ionBackNum;
            totalIonNum += ni;
    }
    //hze1->Write();

    f->Close();
    
    std::ofstream fichier("IBF_summary.txt", std::ios::out | std::ios::app);  //déclaration du flux et ouverture du fichier
        
    if(fichier)  // si l'ouverture a réussi
    {
        fichier << "Vmesh = " << hvMesh << "V " << std::endl;
        fichier << "Total number of ions : " << totalIonNum << std::endl;
        fichier << "Ions that escaped : " << totalIonBackNum << std::endl;
        fichier << "IBF = " << totalIonBackNum * 100./totalIonNum << " %\n\n" << std::endl;
        fichier.close();  // on referme le fichier
    }
    else  { // sinon
        std::cerr << "Erreur à l'ouverture !" << std::endl;
    }
    time_t t1 = time(NULL);
    
    std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
}



