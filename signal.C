#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include "_Utils.C"
#include "parameters.C"

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
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //const int modelNum = 14;
    const bool computeIBF = true;  // if false, it will only compute the number of amplification electrons in the avalanche
    const int nEvents = 1000;  // number of avalanches to simulate
    //____________________
    
    time_t t0 = time(NULL);
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    TString fOutputName;
    
    // Make a gas medium.
    MediumMagboltz* gas = InitiateGas(gasName);
    // Load field map
    ComponentComsol* fm;
    
    int hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvGemDown = 0, hvGemUp = 0, hvMeshTop = 0, hvDrift = 0;
    int saveNum;
    if (modelNum == 1) {
        if (argc != 4) {
            std::cout << "Please enter HVmesh like this: ./signal $hvMesh $hvDrift $saveNum " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDrift = atoi(argv[2]);
        saveNum = atoi(argv[3]);
        fm = InitiateField(modelNum, hvMesh, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDrift, saveNum);
    }
    else if ((modelNum >= 2 && modelNum < 5) || modelNum == 14) {
        if (argc != 5) {
            std::cout << "Please enter HVmesh like this: ./signal $hvDmDown $hvDmUp $hvDrift $saveNum " << std::endl;
            return 0;
        }
        hvDmDown = atoi(argv[1]);
        hvDmUp = atoi(argv[2]);
        hvDrift = atoi(argv[3]);
        saveNum = atoi(argv[4]);
        fm = InitiateField(modelNum, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift, saveNum);
    }
    else if (modelNum >= 5 && modelNum < 8) {
        if (argc != 6) {
            std::cout << "Please enter HVmesh like this: ./signal $hvMesh $hvDmDown $hvDmUp $hvDrift $saveNum " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvDmDown = atoi(argv[2]);
        hvDmUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        saveNum = atoi(argv[5]);
        fm = InitiateField(modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDmDown, hvDmUp, hvDrift, saveNum);
    }
    else if (modelNum >= 8 && modelNum < 10) {
        if (argc != 6) {
            std::cout << "Please enter HVmesh like this: ./signal $hvMesh $hvGemDown $hvGemUp $hvDrift $saveNum " << std::endl;
            return 0;
        }
        hvMesh = atoi(argv[1]);
        hvGemDown = atoi(argv[2]);
        hvGemUp = atoi(argv[3]);
        hvDrift = atoi(argv[4]);
        saveNum = atoi(argv[5]);
        fm = InitiateField(modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift, gas);
        fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvGemDown, hvGemUp, hvDrift, saveNum);
    }
	else if (modelNum >= 10 && modelNum < 12) {
		if (argc != 7) {
			std::cout << "Please enter HVmesh like this: ./signal $hvMesh $hvGemDown $hvGemUp *hvMeshTop $hvDrift $saveNum " << std::endl;
			return 0;
		}
		hvMesh = atoi(argv[1]);
		hvGemDown = atoi(argv[2]);
		hvGemUp = atoi(argv[3]);
		hvMeshTop = atoi(argv[4]);
		hvDrift = atoi(argv[5]);
		saveNum = atoi(argv[6]);
		fm = InitiateField(modelNum, hvMesh, hvGemDown, hvGemUp, hvMeshTop, hvDrift, gas);
		fOutputName = Form("rootFiles/%s/model%d/signal-%d-%d-%d-%d-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvGemDown, hvGemUp, hvMeshTop, hvDrift, saveNum);
	}
    else {std::cout << "Wrong model number" << std::endl; return 0;}
    
    if (!fm || fm->GetMedium(0,0,0) == nullptr) {
        std::cout << "Component COMSOL was not initialized, please fix this" << std::endl;
        return 0;
    }
    
        //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
	int electrodeNum = 0;
	electrodeNum = GetElectrodeNum(modelNum);
	if (electrodeNum == 0) {std::cout << "Warning! Number of electrodes = 0" << std::endl; return 0;}
    //std::cout << damp << " " << width << " " << depth << " " << ddrift << std::endl;
    //return 0;
    
    // Make a sensor.
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
    sensor->SetArea(0, 0, 0, width, depth, ddrift);
    //sensor->SetArea(pitch, pitch, damp-pitch, 3*pitch, 3*pitch, damp+pitch);
    for (int k = 0; k < electrodeNum; k++) sensor->AddElectrode(fm, Form("V%d", k+2));
    //return 0;

    
    // Create ROOT histograms of the signal and a file in which to store them.
    TFile* f = new TFile(fOutputName, "RECREATE");
    //TFile* f = new TFile("rootFiles/test.root", "RECREATE");
    
    TTree *tAvalanche = new TTree("tAvalanche","Gain");
    int ne2 = 0, nWinners = 0;
    tAvalanche->Branch("amplificationElectrons", &nWinners, "amplificationElectrons/I");
	tAvalanche->Branch("avalancheSize", &ne2, "avalancheSize/I");
    Double_t ibfRatio = 0.;
    if (computeIBF) tAvalanche->Branch("ibfRatio", &ibfRatio, "ibfRatio/D");
    
    // Set the signal binning.
    //const int nTimeBins = 10000;
    const double tStep = 0.1;   //ns
    //const double timespace = 1./rate*1.e9;    // in ns
    const double timespace = 1000.;    // in ns // 1µs so that events don't overlap
    const double tStart =  0.;
    //const double tEnd = int(nEvents * timespace + ionDelay);
    const double tEnd = int(nEvents * timespace);
    //const double tStep = (tEnd - tStart) / nTimeBins;
    const int nTimeBins = (tEnd - tStart)/tStep;
    sensor->SetTimeWindow(tStart, tStep, nTimeBins);
    
    // To look at the avalanche
    ViewDrift* driftView = new ViewDrift();
    //driftView->SetArea(390*pitch, 390*pitch, 0, 410*pitch, 410*pitch, damp+6*pitch);
    driftView->SetArea(0, 0, 0, width, depth, ddrift);
    
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
    
    //return 0;
    for (unsigned int i = 0; i < nEvents; ++i) {
        if (i % division == 0) {
            std::cout << "\n\n\n\n" << i << "/" << nEvents << std::endl;
            time_t t = time(NULL);
            PrintTime(t0, t);
            tAvalanche->Write("", TObject::kOverwrite);
        }
        // Initial coordinates of the photon.
        double x0 = width/2. + RndmUniform() * pitch;
        //double y0 = RndmUniform() * depth;
        double y0 = depth/2. + RndmUniform() * pitch;
        double z0 = damp + 2*radius + (ddrift-damp-2*radius)*RndmUniform();
        //double t0 = ( i + RndmUniform() )* timespace;
        double t0 = i * timespace;
        double e = 0;
        aval->AvalancheElectron(x0, y0, z0, t0, e, 0, 0, -1);
        //int ne2 = 0, ni = 0;
		int ni = 0;
        aval->GetAvalancheSize(ne2, ni);
        std::cout << "\nAvalanche size = " << ne2 << std::endl;
        //if (ne2 < 4) {i--; continue;}
        /*
        if (modelNum == 1) { nWinners = ne2;
            //hElectrons->Fill(ne2);    // ok for modelNum < 4
            //continue;
        }
         */
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        double xi1, yi1, zi1, ti1;
        double xi2, yi2, zi2, ti2;
        int status;
        int ionBackNum = 0;
        nWinners = 0;
        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status);
            if (ze2 < 0.01) nWinners++;
            if (!computeIBF) continue;
            drift->DriftIon(xe1, ye1, ze1, te1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
            if (zi2 > 1.2*damp) ionBackNum+=1;
        }
        std::cout << "nWinners = " << nWinners << " / " << ne2 << std::endl;
		if (computeIBF) ibfRatio = (double)ionBackNum/ni;
            //std::cout << ibfRatio << std::endl;
		tAvalanche->Fill();
        //hTransparencySA->Fill(electronsBelowSA*1./electronsAboveSA);
    }
    tAvalanche->Write("", TObject::kOverwrite);
    
    //std::cout << "induced charge for electrode " << 2 << " = " << sensor->GetInducedCharge("V2") << std::endl;

    if (computeIBF) {
        for (int k = 0; k < electrodeNum; k++) {
            TTree *tSignal = new TTree(Form("tSignal_%d",k+2),"Currents");
            Double_t ft = 0., fct = 0., fce = 0., fci = 0.;
            Double_t fctInt = 0, fceInt = 0, fciInt = 0;
            tSignal->Branch("time", &ft, "time/D");
            tSignal->Branch("totalCurrent", &fct, "totalCurrent/D");
            tSignal->Branch("electronCurrent", &fce, "electronCurrent/D");
            tSignal->Branch("ionCurrent", &fci, "ionCurrent/D");
            tSignal->Branch("totalCurrentInt",&fctInt,"totalCurrentInt/D");     // Integral
            tSignal->Branch("electronCurrentInt",&fceInt,"electronCurrentInt/D");
            tSignal->Branch("ionCurrentInt",&fciInt,"ionCurrentInt/D");
            
            TTree *tCharge = new TTree(Form("tInducedCharge_%d",k+2),"InducedCharge");
            Double_t tic = 0., eic = 0., iic = 0.;
            tCharge->Branch("totalInducedCharge", &tic, "totalInducedCharge/D");
            tCharge->Branch("electronInducedCharge", &eic, "electronInducedCharge/D");
            tCharge->Branch("ionInducedCharge", &iic, "ionInducedCharge/D");
            
     
            //int evNum = 0;
            for (int j = 0; j < nTimeBins; j++) {
                ft = j * tStep;
                //std::cout << ft << std::endl;
                fct = sensor->GetSignal(Form("V%d", k+2), j) / ElementaryCharge;
                fce = sensor->GetElectronSignal(Form("V%d", k+2), j) / ElementaryCharge;
                fci = sensor->GetIonSignal(Form("V%d", k+2), j) / ElementaryCharge;
                fctInt+=fct*tStep;
                fceInt+=fce*tStep;
                fciInt+=fci*tStep;
                tSignal->Fill();
                
                tic+=fct*tStep;
                eic+=fce*tStep;
                iic+=fci*tStep;
                
                if (int((j+1) * tStep/timespace) == (j+1) * tStep/timespace) {
                    //tAvalanche->GetEntry(evNum);
                    //std::cout << "number of amplification electrons = " << nWinners << std::endl;
                    //std::cout << "total charge induced = " << int(tic) << std::endl;
                    //evNum++;
                    tCharge->Fill();
                    tic = 0; eic = 0; iic = 0;
                }
            }
            tSignal->Write("", TObject::kOverwrite);
            tCharge->Write("", TObject::kOverwrite);
            //std::cout << "induced charge for electrode " << k+2 << " = " << sensor->GetInducedCharge(Form("V%d", k+2)) << std::endl;
        }
    }
    f->Close();

    
    // Plot the signal.
    const bool plotSignal = false;
    if (plotSignal) {
        TCanvas* cSignal = new TCanvas("signal", "Signal");
        ViewSignal* vSignal = new ViewSignal();
        gPad->SetLeftMargin(0.15);
        vSignal->SetSensor(sensor);
        vSignal->SetCanvas(cSignal);
        //void PlotSignal(const std::string& label, const bool total = true, const bool electron = false, const bool ion = false);
        vSignal->SetRangeX(0, nEvents*timespace+1000);
        vSignal->PlotSignal("V2", true, true, true);    // signal de la mesh
        cSignal->SaveAs("Figures/Signal.pdf");
    }
    
    
    time_t t1 = time(NULL);
    
    std::cout << "\n" << nEvents << " events simulated" << std::endl;
    PrintTime(t0, t1);
    
    //app.Run(true);
}


