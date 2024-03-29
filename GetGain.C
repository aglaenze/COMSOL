#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>

#include "Include/Data.C"
#include "Include/Spectrum.C"

/*
 TO-DO:
 - draw gain as a function of field ratio
 */


int GetGain(bool ibf = true) {
    
    //______________________
    // variables
    std::string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    //std::string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 23;
    //____________________
    
    time_t t0 = time(NULL);
    gStyle->SetTitleFontSize(.06);
    gStyle->SetTitleSize(.06);
    
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetLabelSize(.05, "XY");
    //gStyle->SetMarkerSize(0.3);
    gStyle->SetMarkerStyle(20);
    gStyle->SetTextSize(.05);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    
    const TString path = Form("rootFiles/%s/model%d/", gasName.c_str(), modelNum);
    
    cout << path << endl;
    
    double damp = 0.0128, ddrift = 0.5;
    //if (modelNum == 22) damp = 0.0125;
    
    std::map <std::string, int, NoSorting> electrode;
    LoadElectrodeMap(modelNum, electrode);
    
    int readoutElectrode = 0;
    int driftElectrode = 0;
    
    std::map<std::string, int>::iterator it = electrode.begin();
    for (it=electrode.begin(); it!=electrode.end(); ++it) {
        //std::cout << it->first << " => " << it->second << '\n';
        if (it->first == "pad") readoutElectrode = it->second;    // could be mesh, it depends on where you want to read
        else if (it->first == "drift") driftElectrode = it->second;
    }
    if (readoutElectrode == 0 || driftElectrode == 0) {std::cout << "Did not find drift or pad electrode" << std::endl; return 0;}
    
    
    // Get number of files to look at
    //Int_t num = 1;
    string filesName = "fe-spectrum-convoluted";
    if (!ibf) filesName += "-noibf";
    Int_t num = GetNumberOfFiles(path, filesName);
    //num -= 1;
    Int_t num2 = int(num/2.);
    if ((double)num2 < (double)num/2.) num2+=1;
    
    TCanvas* c2 = new TCanvas("c2","c2", 600, 300*num2);
    //c2->Divide(2);
    c2->Divide(2, num2);
    
    //num = 4;
    Double_t hvMeshList[num], hvDriftList[num], hvRatioList[num];
    Double_t fieldList[num];
    Double_t gainList[num], gainFeChargeList[num], gainErrorList[num], gainFeChargeErrorList[num];
    for (unsigned int k = 0; k < num; k++) {
        Int_t hvMesh = 0, hvDmDown = 0, hvDmUp = 0, hvGemDown = 0, hvGemUp = 0, hvDrift = 0;
        TString fileName;
        //if (modelNum == 1 || modelNum == 16 || modelNum == 17 || modelNum == 18 || modelNum == 22 || modelNum == 23) {
        int step = 20;
        hvMesh = 340+step*k;
        hvDrift = 540+step*k;
        fileName = Form("rootFiles/%s/model%d/fe-spectrum-convoluted-%d-%d.root", gasName.c_str(), modelNum, hvMesh, hvDrift);
        //}
        /*
         else if (modelNum >= 2) {
         hvDmDown = 300;
         hvDmUp = 600;
         hvDrift = 800;
         fileName = Form("rootFiles/%s/model%d/fe-spectrum-convoluted-%d-%d-%d.root", gasName.c_str(), modelNum, hvDmDown, hvDmUp, hvDrift);
         }
         */
        hvMeshList[k] = hvMesh;
        hvDriftList[k] = hvDrift;
        hvRatioList[k] = (double)hvMesh/(hvDrift-hvMesh)*(ddrift-damp)/damp;
        //hvRatioList[k] = (double)hvMesh/(hvDrift-hvMesh)*ddrift/damp;
        fieldList[k] = (double)hvMesh/damp/1000.;
        
        c2->cd(k+1);
        Double_t gain = 0, gainError = 0;
        DrawFeConvolution(fileName, gain, gainError);
        gainList[k] = gain;
        gainErrorList[k] = gainError;
        
    }
    c2->SaveAs(Form("Figures/model%d/gains-%s.pdf", modelNum, gasName.c_str()));
    
    // Load data
    const Int_t dataNum = dataQuantity(gasName);
    Double_t hvMeshListData[dataNum], hvDriftListData[dataNum], gainListData[dataNum], gainErrorListData_old[dataNum], gainListData_old[dataNum], gainErrorListData[dataNum];
    Double_t fieldListData[dataNum];
    LoadGainData(gasName, dataNum, hvMeshListData, hvDriftListData, gainListData, gainErrorListData, gainListData_old, gainErrorListData_old);
    Double_t hvRatioListData[dataNum];
    for (int i = 0; i < dataNum; i++) {
        hvRatioListData[i] = hvMeshListData[i]/(hvDriftListData[i]-hvMeshListData[i]) * (ddrift-damp)/damp;
        fieldListData[i] = hvMeshListData[i]/damp;
    }
    
    // Draw results of simulation
    TCanvas* c3 = new TCanvas("c3", "c3", 600, 300);
    c3->Divide(2);
    
    c3->cd(1);
    gPad->SetLogy();
    gPad->SetGrid();
    TGraphErrors* grSim1 = new TGraphErrors(num, fieldList, gainList, 0, gainErrorList);
    TGraphErrors* grSim2 = new TGraphErrors(num, fieldList, gainFeChargeList, 0, gainFeChargeErrorList);
    grSim1->SetTitle("Gain curve in the Micromegas");
    //grSim1->GetXaxis()->SetTitle( "V_{mesh}" );
    grSim1->GetXaxis()->SetTitle( "E_{amp} (kV/cm)" );
    grSim1->GetYaxis()->SetTitle( "Gain" );
    //grSim1->GetHistogram()->GetXaxis()->SetLimits(330, 450);
    //grSim1->GetHistogram()->GetXaxis()->SetLimits(20000, 40000);
    grSim1->GetHistogram()->SetMinimum(2.e2);   // along Y axis
    grSim1->GetHistogram()->SetMaximum(8.e4);   // along Y axis
    grSim1->SetMarkerStyle(20);
    grSim1->Draw("AP");
    grSim2->SetMarkerStyle(20);
    //grSim2->Draw("P same");
    
    // Fit the gain curve with an exponential function
    TF1* fExp = new TF1( "FitFunctionExp", FitFunctionExp, fieldList[0]-5, fieldList[num-1]+5, 2 );
    fExp->SetParName( 0, "const" );
    fExp->SetParName( 1, "slope" );
    grSim1->Fit( fExp, "0");
    fExp->SetLineColor(2);
    fExp->Draw("same");
    PutText( 0.2, 0.7, Form( "#bf{y_{simulation} = exp( %.3f + %.3f x)}", fExp->GetParameter(0), fExp->GetParameter(1) ) );
    
    // Same with data
    
    TGraphErrors* grData1 = new TGraphErrors(dataNum, fieldListData, gainListData, 0, gainErrorListData);
    //grData1->Draw("P same");
    TGraphErrors* grData2 = new TGraphErrors(dataNum, fieldListData, gainListData_old, 0, gainErrorListData_old);
    //grData2->Draw("P same");
    
    TF1* fExp2 = new TF1( "FitFunctionExp2", FitFunctionExp, hvMeshListData[0]-5, hvMeshListData[dataNum-1]+5, 2 );
    fExp2->SetParName( 0, "const" );
    fExp2->SetParName( 1, "slope" );
    fExp2->SetLineColor(3);
    grData1->Fit( fExp2, "0" );
    //fExp2->Draw("same");
    //PutText( 0.2, 0.8, Form( "#bf{y_{data} = exp( %.3g + %.3g x)}", fExp2->GetParameter(0), fExp2->GetParameter(1) ) );
    
    TF1* fExp3 = new TF1( "FitFunctionExp3", FitFunctionExp, hvMeshListData[0]-5, hvMeshListData[dataNum-1]+5, 2 );
    fExp3->SetParName( 0, "const" );
    fExp3->SetParName( 1, "slope" );
    fExp3->SetLineColor(4);
    grData2->Fit( fExp3, "0" );
    //fExp3->Draw("same");
    //PutText( 0.2, 0.75, Form( "#bf{y_{data old} = exp( %.3g + %.3g x)}", fExp3->GetParameter(0), fExp3->GetParameter(1) ) );
    //PutText( 0.2, 0.75, Form( "#bf{y_{data} = exp( %.3g + %.3g x)}", fExp3->GetParameter(0), fExp3->GetParameter(1) ) );
    
    TLegend* legend = new TLegend(0.7,0.2,0.9,0.4);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    /*
     legend->AddEntry(fExp2,"Data", "l");
     legend->AddEntry(fExp3,"Data (old)", "l");
     */
    legend->AddEntry(fExp3,"Data", "l");
    legend->AddEntry(fExp,"Simulation", "l");
    //legend->Draw();
    
    
    c3->cd(2);
    gPad->SetLogy();
    gPad->SetGrid();
    TGraphErrors* grSimRatio1 = new TGraphErrors(num, hvRatioList, gainList, 0, gainErrorList);
    TGraphErrors* grSimRatio2 = new TGraphErrors(num, hvRatioList, gainFeChargeList, 0, gainFeChargeErrorList);
    grSimRatio1->SetTitle("Gain curve in the Micromegas");
    grSimRatio1->GetXaxis()->SetTitle( "E_{amp}/E_{drift}" );
    grSimRatio1->GetYaxis()->SetTitle( "Gain" );
    grSimRatio1->GetHistogram()->GetXaxis()->SetLimits(hvRatioList[0]-5, hvRatioList[num-1]+5);
    grSimRatio1->GetHistogram()->SetMinimum(2.e2);   // along Y axis
    grSimRatio1->GetHistogram()->SetMaximum(8.e4);   // along Y axis
    grSimRatio1->SetMarkerStyle(20);
    grSimRatio1->Draw("AP");
    grSimRatio2->SetMarkerStyle(20);
    //grSimRatio2->Draw("P same");
    
    // Fit the gain curve with an exponential function
    TF1* fExpRatio = new TF1( "FitFunctionExp", FitFunctionExp, hvRatioList[0]-2, hvRatioList[num-1]+2, 2 );
    fExpRatio->SetParName( 0, "const" );
    fExpRatio->SetParName( 1, "slope" );
    grSimRatio1->Fit( fExpRatio, "0");
    fExpRatio->SetLineColor(2);
    fExpRatio->Draw("same");
    PutText( 0.2, 0.7, Form( "#bf{y_{simulation} = exp( %.3f + %.3f x)}", fExpRatio->GetParameter(0), fExpRatio->GetParameter(1) ) );
    
    TGraphErrors* grDataRatio1 = new TGraphErrors(dataNum, hvRatioListData, gainListData, 0, gainErrorListData);
    grDataRatio1->SetMarkerStyle(20);
    //grDataRatio1->Draw("P same");
    TGraphErrors* grDataRatio2 = new TGraphErrors(dataNum, hvRatioListData, gainListData_old, 0, gainErrorListData_old);
    grDataRatio2->SetMarkerStyle(20);
    //grDataRatio2->Draw("P same");
    
    TF1* fExpRatio2 = new TF1( "FitFunctionExp2", FitFunctionExp, hvRatioListData[0]-2, hvRatioListData[dataNum-1]+2, 2 );
    fExpRatio2->SetParName( 0, "const" );
    fExpRatio2->SetParName( 1, "slope" );
    fExpRatio2->SetLineColor(3);
    grDataRatio1->Fit( fExpRatio2, "0" );
    //fExpRatio2->Draw("same");
    //PutText( 0.2, 0.8, Form( "#bf{y_{data} = exp( %.3g + %.3g x)}", fExpRatio2->GetParameter(0), fExpRatio2->GetParameter(1) ) );
    
    TF1* fExpRatio3 = new TF1( "FitFunctionExp3", FitFunctionExp, hvRatioListData[0]-2, hvRatioListData[dataNum-1]+2, 2 );
    fExpRatio3->SetParName( 0, "const" );
    fExpRatio3->SetParName( 1, "slope" );
    fExpRatio3->SetLineColor(4);
    grDataRatio2->Fit( fExpRatio3, "0" );
    //fExpRatio3->Draw("same");
    //PutText( 0.2, 0.75, Form( "#bf{y_{data old} = exp( %.3g + %.3g x)}", fExpRatio3->GetParameter(0), fExpRatio3->GetParameter(1) ) );
    //PutText( 0.2, 0.75, Form( "#bf{y_{data} = exp( %.3g + %.3g x)}", fExpRatio3->GetParameter(0), fExpRatio3->GetParameter(1) ) );
    
    TLegend* legendRatio = new TLegend(0.7,0.2,0.9,0.4);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    /*
     legendRatio->AddEntry(fExpRatio2,"Data", "l");
     legendRatio->AddEntry(fExpRatio3,"Data (old)", "l");
     */
    
    //legendRatio->AddEntry(fExpRatio3,"Data", "l");
    legendRatio->AddEntry(fExpRatio,"Simulation", "l");
    //legendRatio->Draw();
    
    c3->SaveAs(Form("Figures/model%d/GainCurve-%s.pdf", modelNum, gasName.c_str()));
    
    
    time_t t1 = time(NULL);
    //PrintTime(t0, t1);
    
    // Prints lists
    vector<double> gainVec = {};
    vector<double> gainErrorVec = {};
    vector<double> fieldVec = {};
    vector<double> fieldRatioVec = {};
    for (int k = 0; k < num; k++) {
        gainVec.push_back(gainList[k]);
        gainErrorVec.push_back(gainErrorList[k]);
        fieldVec.push_back(fieldList[k]);
        fieldRatioVec.push_back(hvRatioList[k]);
    }
    PrintListVec( "gainListSim", gainVec, "%.3f" );
    PrintListVec( "gainErrorListSim", gainErrorVec, "%.3f" );
    PrintListVec( "fieldListSim", fieldVec, "%.2f" );
    PrintListVec( "fieldRatioListSim", fieldRatioVec, "%.2f" );
    
    return 0;
}


