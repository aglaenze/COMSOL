#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>

#include "Include/Ibf.C"
#include "Include/Electrons.C"


int Analyse(int modelNum, std::string gasName, std::vector<int> hvList) {
    
    bool ionAnalysis = false;
    
    time_t t0 = time(NULL);
    
    LoadStyle();
    
    int electrodeNum = GetElectrodeNum(modelNum);
    if ((int)hvList.size() != electrodeNum-1) {std::cout << "Wrong hv input" << std::endl; return 0;}
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    // input and output files
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    TString fileName = path + "fe-spectrum-convoluted";
    TString fSignalName = path + "signal";
    TString outputName = Form("Figures/model%d/analysis-%s", modelNum, gasName.c_str());
    TString outputName2 = Form("Figures/model%d/ion-analysis-%s", modelNum, gasName.c_str());
    for (int k = 0; k< (int)hvList.size(); k++) {
        fileName += Form("-%d", hvList[k]);
        fSignalName += Form("-%d", hvList[k]);
        outputName += Form("-%d", hvList[k]);
        outputName2 += Form("-%d", hvList[k]);
    }
    fileName += ".root";
    fSignalName += ".root";
    outputName+=".pdf";
    outputName2+=".pdf";
    //outputName+=".root";
    
    std::cout << fSignalName << std::endl;
    
    
    
    
    /*
     std::vector <Int_t> hv = {};
     TObjArray* matches = TPRegexp( "-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?-?(\\d+)?\\.root$" ).MatchS( fileName );
     for (int k = 0; k<matches->GetLast(); k++) {
     hv.push_back( (static_cast<TObjString*>(matches->At(k+1)))->GetString().Atoi() );
     std::cout << hv[k] << std::endl;
     }
     */
    
    TCanvas* cv = new TCanvas("cv","cv", 1200, 1000);
    //cv->Divide(2);
    cv->Divide(3, 2);
    
    // Start with drawing the transparency
    cv->cd(4);
    DrawTransparency(modelNum, fSignalName);
    //DrawGains(modelNum, fileName);
    //return 0;
    
    // now draw the gain with the 2 different ways
    
    TFile* fSignal = TFile::Open(fSignalName, "READ");
    TTree* tAvalanche = (TTree*) fSignal->Get("tAvalanche");
    
    
    // Now draw both spectra
    cv->cd(2);
    
    //DrawAmplificationElectrons(gasName, fSignalName, false);
    // Ignore convolution in the end
    DrawFeConvolution(fileName);
    string readout = "mesh";
    //DrawFeChargeConvolution(modelNum, fileName, readout);
    
    
    std::cout << "\n\nStarting to draw the IBF now\n\n" << std::endl;
    
    // And finally draw the IBF
    
    cv->cd(3);
    //DrawIbf(modelNum, fSignalName);
    DrawConvolutedIbf(fileName);
    
    
    cv->cd(1);
    DrawDetector(modelNum, hvList);
    
    TText* txtGas = new TText(.4,.95,Form("Gas: %s", gasName.c_str()));
    txtGas->Draw("same");
    
    // Draw distribution of where electrons are created
    cv->cd(5);
    gPad->SetLogy();
    double zmin = 0;
    double zmax = damp*1.2;
    //DrawElectronsStart(modelNum, fSignalName, zmin, zmax);
    //DrawElectrodes(modelNum, zmin, zmax, true);
    
    DrawElectrodeStart(modelNum, fSignalName);
    
    cv->cd(6);
    gPad->SetLogy();
    gPad->SetBottomMargin(0.1);
    //DrawDyingIons(modelNum, fSignalName);
    DrawDyingIonsOverlaid(modelNum, fSignalName);
    
    cv->SaveAs(outputName);
    
    if (!ionAnalysis) {
        time_t t1 = time(NULL);
        PrintTime(t0, t1);
        
        return 0;
    }
    
    // Draw histograms of IBF (see in Transparency.C)
    
    TCanvas* cv2 = new TCanvas("cv2","cv2", 1200, 1000);
    //cv->Divide(2);
    cv2->Divide(3, 2);
    
    cv2->cd(1);
    DrawNionsInDriftRegion(fSignalName);
    
    cv2->cd(2);
    DrawDyingIons(modelNum, fSignalName);
    cv2->cd(3);
    DrawDyingIons2d(modelNum, fSignalName);
    cv2->cd(4);
    DrawDyingIons3d(modelNum, fSignalName);
    
    cv2->SaveAs(outputName2);
    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    return 0;
}


