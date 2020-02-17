#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>

#include "_Utils.C"
#include "parameters.C"

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
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 10;
    //____________________
    
    time_t t0 = time(NULL);
    
    char* fOutputName;
    
    // Make a gas medium.
    MediumMagboltz* gas = InitiateGas(gasName);
    // Load field map
    ComponentComsol* fm;
    
    int hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvDrift = 0;
    if (modelNum == 1) {
        if (argc != 3 ) {
            std::cout << "Please enter HVmesh like this: ./gain $hvMesh" << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDrift = atoi(argv[2]);
        fm = InitiateField(modelNum, hvMesh, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/ibf-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDrift);
    }
    else if (modelNum > 6 && modelNum < 10) {
        if (argc != 4) {
            std::cout << "Please enter HVmesh like this: ./gain $hvDmDown $hvDmUp" << std::endl;
            return 0;
        }
        hvDmDown = atoi(argv[1]);
        hvDmUp = atoi(argv[2]);
        hvDrift = atoi(argv[3]);
        fm = InitiateField(modelNum, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/ibf-%d-%d-%d.root", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift);
    }
    else if (modelNum == 10) {
        if (argc != 5) {std::cout << "Please enter HVmesh like this: ./gain $hvMesh $hvDmDown $hvDmUp " << std::endl; return 0;}
        hvMesh = atoi(argv[1]);
        hvDmDown = atoi(argv[2]);
        hvDmUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        fm = InitiateField(modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/ibf-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift);
    }
    else {std::cout << "Wrong model number" << std::endl; return 0;}
    
    std::cout << fOutputName << std::endl;
    //return 0;

    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);

    /*
    Medium* wire = new Medium();
    wire->SetTemperature(293.15);
    wire->SetPressure(AtmosphericPressure);
    //wire->SetComposition("Stainless steel");
     */
    
    
    // Make a sensor.
    Sensor sensor;
    sensor.AddComponent(fm);
    sensor.SetArea(0, 0, 0, width, depth, ddrift);
    //sensor.SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
    
    return 0;

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
    TH1F* hIbf = new TH1F("hIbf", "hIbf", 10000, 0, 100);
    //TH1::StatOverflows(true);
    hIbf->SetXTitle("IBF (%)");
    hIbf->SetYTitle("# counts");
    /*
    TH1F* hTransparencySA = new TH1F("hTransparencySA", "hTransparencySA", 1000, 0., 1.);
    hTransparencySA->SetYTitle("# counts");
    hTransparencySA->SetXTitle("fraction of e- that passed");
     */
    
    //TH1F* hze1 = new TH1F("hze1", "hze1", 10000, 0, damp);

    // Write the histograms to the TFile.
    TFile* f = new TFile(fOutputName, "RECREATE");
    
    const int nEvents = 1000;
    int division = int(nEvents/20);
    
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (i % 100 == 0) {
            std::cout << "\n\n\n\n" << i << "/" << nEvents << "\n";
            time_t t = time(NULL);
            PrintTime(t0, t);
        }
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
        // hElectrons->Fill(ne2);
        // to look at avalanche
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        double xi1, yi1, zi1, ti1;
        double xi2, yi2, zi2, ti2;
        int status;
        int ionBackNum = 0;
        //int electronsAboveSA = 0;
        //int electronsBelowSA = 0;
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
            /*
            if (ze1 > 0.0147) {
                electronsAboveSA++;
                if (ze2 < 0.008) electronsBelowSA++;
            }
             */
            //std::cout << "departure of the electron in x y z : " << xe1 << " " << ye1 << " " <<  ze1 << std::endl;
            //std::cout << "arrival of the electron in x y z : " << xe2 << " " << ye2 << " " <<  ze2 << std::endl;
            drift->DriftIon(xe1, ye1, ze1, te1);
            //hze1->Fill(ze1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
            //std::cout << "arrival of the ion in z : " << zi2 << std::endl;
            if (zi2 > 3*damp) ionBackNum+=1;
        }
        hIbf->Fill( ionBackNum*100./ni);
        //hTransparencySA->Fill(electronsBelowSA*1./electronsAboveSA);
        if (i % division == 0) {
            hIbf->Write("", TObject::kOverwrite);
            //hTransparencySA->Write("", TObject::kOverwrite);
        }
        totalIonBackNum += ionBackNum;
        totalIonNum += ni;
    }
    //hze1->Write();

    f->Close();
    
    std::ofstream fichier("IBF_summary.txt", std::ios::out | std::ios::app);  //déclaration du flux et ouverture du fichier
        
    if(fichier)  // si l'ouverture a réussi
    {
        //fichier << "Vmesh = " << hvMesh << "V " << std::endl;
        fichier << "VmeshDown = " << hvDmDown << "V " << std::endl;
        fichier << "VmeshUp = " << hvDmUp << "V " << std::endl;
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



