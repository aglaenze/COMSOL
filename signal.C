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
#include "Garfield/ViewSignal.hh"
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
        std::cout << "Please enter HVmesh like this: ./signal $hvMesh " << std::endl;
        return 0;
    }
    const int hvMesh = atoi(argv[1]);
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    // Set up detector geometry
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
    
    // Load the field map.
    std::string dataFolder = Form("COMSOL_data/model%d/", modelNum);
    std::string dataFile = dataFolder + Form("ewfield_%dV.txt", hvMesh);
    ComponentComsol* fm = new ComponentComsol();
    fm->Initialise(dataFolder+"mesh.mphtxt", dataFolder+"dielectrics.dat", dataFile);
    //fm->SetWeightingField("COMSOL_data/base/wfieldmesh.txt", "V2");
    std::cout << "\n\nWP = " << fm->WeightingPotential(pitch/2, pitch/2., damp+pitch, "V2") << std::endl;
    
    fm->PrintMaterials();
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
    fm->PrintRange();
    
    fm->SetMedium(0, gas);
    
    fm->PrintMaterials();
    
    // Make a sensor.
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
    sensor->SetArea(0, 0, 0, width, depth, ddrift);
    //sensor->SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
    sensor->AddElectrode(fm, "V2");
    sensor->AddElectrode(fm, "V3");
    
    // Set the signal binning.
    const int nEvents = 30;
    const double timespace = 10;    // 1 event every 10 ns
    const double ionDelay = 1.e6;    // time to collect all ions at the drift ˜1ms
    double tEnd = nEvents * timespace + ionDelay;
    int nsBins = int(tEnd);
    sensor->SetTimeWindow(0., tEnd / nsBins, nsBins);
    
    // To look at the avalanche
    ViewDrift* driftView = new ViewDrift();
    driftView->SetArea(390*pitch, 390*pitch, 0, 410*pitch, 410*pitch, damp+6*pitch);
    
    // Create an avalanche object
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    aval->EnablePlotting(driftView);
    aval->EnableSignalCalculation();
    
    AvalancheMC* drift = new AvalancheMC();
    drift->SetSensor(sensor);
    drift->SetDistanceSteps(2.e-4);
    drift->EnablePlotting(driftView);
    drift->EnableSignalCalculation();
    
    int division = int(nEvents/10);
    
    
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (i % division == 0) std::cout << i << "/" << nEvents << "\n";
        // Initial coordinates of the photon.
        double x0 = width/2. + RndmUniform() * pitch;
        //double y0 = RndmUniform() * depth;
        double y0 = depth/2. + RndmUniform() * pitch;
        double z0 = damp + radius + (ddrift-damp-radius)*RndmUniform();
        double t0 = i*timespace;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        int ne2 = 0, ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        std::cout << "\nAvalanche size = " << ne2 << std::endl;
        if (ne2 < 2) {i--; continue;}
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        double xi1, yi1, zi1, ti1;
        double xi2, yi2, zi2, ti2;
        int status;
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
            drift->DriftIon(xe1, ye1, ze1, te1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
        }
    }
    
    const bool extract = true;
    if (extract) {
        // Extract the calculated signal.
        double bscale = tEnd / nsBins;  // time per bin
        double sum = 0.;  // to keep a running sum of the integrated signal
        double sumE = 0.;
        double sumI = 0.;
        double sumd = 0.;  // drift electrode
        double sumEd = 0.;
        double sumId = 0.;
        
        // Create ROOT histograms of the signal and a file in which to store them.
        const char* name = Form("rootFiles/%s/model%d/signal_%dV.root", gasName.c_str(), modelNum, hvMesh);
        TFile* f = new TFile(name, "RECREATE");
        TH1F* hS = new TH1F("hhMesh", "Total signal (mesh)", nsBins, 0, tEnd);
        TH1F* hInt = new TH1F("hIntMesh", "Integrated signal (mesh)", nsBins, 0, tEnd);
        TH1F* hSE = new TH1F("hhEMesh", "Electron signal (mesh)", nsBins, 0, tEnd);
        TH1F* hIntE = new TH1F("hIntEMesh", "Integrated electron signal (mesh)", nsBins, 0, tEnd);
        TH1F* hSI = new TH1F("hhIMesh", "Ion signal (mesh)", nsBins, 0, tEnd);
        TH1F* hIntI = new TH1F("hIntIMesh", "Integrated ion signal (mesh)", nsBins, 0, tEnd);

        TH1F* hSd = new TH1F("hhDrift", "Total signal (drift)", nsBins, 0, tEnd);
        TH1F* hIntd = new TH1F("hIntDrift", "Integrated signal (drift)", nsBins, 0, tEnd);
        TH1F* hSEd = new TH1F("hhEDrift", "Electron signal (drift)", nsBins, 0, tEnd);
        TH1F* hIntEd = new TH1F("hIntEDrift", "Integrated electron signal (drift)", nsBins, 0, tEnd);
        TH1F* hSId = new TH1F("hhIDrift", "Ion signal (drift)", nsBins, 0, tEnd);
        TH1F* hIntId = new TH1F("hIntIDrift", "Integrated ion signal (drift)", nsBins, 0, tEnd);

        // Fill the histograms with the signals.
        //  Note that the signals will be in C/(ns*binWidth), and we will divide by e
        // to give a signal in e/(ns*binWidth).
        //  The total signal is then the integral over all bins multiplied by the bin
        // width in ns.
        for (int i = 0; i < nsBins; i++) {
            double wt = sensor->GetSignal("V2", i) / ElementaryCharge;
            sum += wt;
            hS->Fill(i * bscale, wt);
            hInt->Fill(i * bscale, sum);
            double wtE = sensor->GetElectronSignal("V2", i) / ElementaryCharge;
            sumE += wtE;
            hSE->Fill(i * bscale, wtE);
            hIntE->Fill(i * bscale, sumE);
            double wtI = sensor->GetIonSignal("V2", i) / ElementaryCharge;
            sumI += wtI;
            hSI->Fill(i * bscale, wtI);
            hIntI->Fill(i * bscale, sumI);
            
            double wtd = sensor->GetSignal("V3", i) / ElementaryCharge;
            sumd += wtd;
            hSd->Fill(i * bscale, wtd);
            hIntd->Fill(i * bscale, sumd);
            double wtEd = sensor->GetElectronSignal("V3", i) / ElementaryCharge;
            sumEd += wtEd;
            hSEd->Fill(i * bscale, wtEd);
            hIntEd->Fill(i * bscale, sumEd);
            double wtId = sensor->GetIonSignal("V3", i) / ElementaryCharge;
            sumId += wtId;
            hSId->Fill(i * bscale, wtId);
            hIntId->Fill(i * bscale, sumId);
        }
        //std::cout << "IBF = " << sumd*100./sum - 100./gain << "%" << std::endl;
        std::cout << "IBF = " << sumd*100./sum << "%" << std::endl;
        
        hS->SetOption("hist");
        hInt->SetOption("hist");
        hSE->SetOption("hist");
        hIntE->SetOption("hist");
        hSI->SetOption("hist");
        hIntI->SetOption("hist");
        hSd->SetOption("hist");   // drift electrode
        hIntd->SetOption("hist");
        hSEd->SetOption("hist");
        hIntEd->SetOption("hist");
        hSId->SetOption("hist");
        hIntId->SetOption("hist");
        
        hS->SetXTitle("time (ns)");
        hS->SetYTitle("# charges / unit time");
        hSE->SetXTitle("time (ns)");
        hSE->SetYTitle("# charges / unit time");
        hSI->SetXTitle("time (ns)");
        hSI->SetYTitle("# charges / unit time");
        hSd->SetXTitle("time (ns)");
        hSd->SetYTitle("# charges / unit time");
        hSEd->SetXTitle("time (ns)");
        hSEd->SetYTitle("# charges / unit time");
        hSId->SetXTitle("time (ns)");
        hSId->SetYTitle("# charges / unit time");
        
        
        // Write the histograms to the TFile.
        hS->Write();    // mesh electrode
        hInt->Write();
        hSE->Write();
        hIntE->Write();
        hSI->Write();
        hIntI->Write();

        hSd->Write();   // drift electrode
        hIntd->Write();
        hSEd->Write();
        hIntEd->Write();
        hSId->Write();
        hIntId->Write();

        f->Close();
    }
    
    // Plot the signal.
    const bool plotSignal = true;
    if (plotSignal) {
        TCanvas* cSignal = new TCanvas("signal", "Signal");
        ViewSignal* vSignal = new ViewSignal();
        vSignal->SetSensor(sensor);
        vSignal->SetCanvas(cSignal);
        //void PlotSignal(const std::string& label, const bool total = true, const bool electron = false, const bool ion = false);
        vSignal->SetRangeX(0, nEvents*timespace+200);
        vSignal->PlotSignal("V2", true, true, true);
        cSignal->SaveAs("Figures/Signal.pdf");
    }
    
    
    time_t t1 = time(NULL);
    
    std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
}


