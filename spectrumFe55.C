/*
Note that the spectrum shape only depends on the detector geometry
 Photons create electrons with energy = 0
*/

#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "Include/Utils.C"
#include "Include/Geometry.C"
#include "Include/Initiate.C"

#include "Garfield/ComponentComsol.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/TrackHeed.hh"
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
using namespace std;

int main(int argc, char * argv[]) {
    
    //______________________
    // variables
    //string gasName = "Ar-CO2"; // Ar-iC4H10 or Ne or Ar-CO2
    string gasName = "Ar-iC4H10"; // Ar-iC4H10 or Ne or Ar-CO2
    const int modelNum = 1;
    const bool variableThickness = false;   // to see how many photons converted depending on the thickness of the detector
    //____________________
    
    time_t t0 = time(NULL);
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    MediumMagboltz* gas = InitiateGas(gasName);
    
    //Load geometry parameters
    double damp = 0., ddrift = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    LoadParameters(modelNum, damp, ddrift, radius, pitch, width, depth);
    
    
    if (variableThickness) {
        const char* nameOut = Form("rootFiles/%s/ConversionNumber.root", gasName.c_str());
        TFile* fOut = new TFile(nameOut, "RECREATE");
        
        double ddriftMin = 0.3;
        double zStep = 0.2;
        int nStep = 50;
        TH1F* hConversion = new TH1F("hConversion", "Proportion of photons that converted", nStep, ddriftMin-zStep/2, ddriftMin+nStep*zStep-zStep/2);

        for (int k = 0; k<nStep; k++) {
            ddrift = ddriftMin+k*zStep; // cm

            cout << "\n\n\n" << nStep-k << " to go... " << endl;
            time_t t1 = time(NULL);
            PrintTime(t0, t1);
            
            // Set up detector geometry
            GeometrySimple* geo = new GeometrySimple();
            // Set up the bulk
            SolidBox bulk = SolidBox(width/2., depth/2., ddrift/2., width/2., depth/2., ddrift/2.);
            geo->AddSolid(&bulk, gas);
            
            // Make a component
            ComponentConstant field;
            field.SetGeometry(geo);

            // Make a sensor.
            Sensor sensor;
            sensor.AddComponent(&field);

            // Use Heed for simulating the photon absorption.
            TrackHeed track;
            track.SetSensor(&sensor);
            //track.EnableElectricField();
            // Histogram
            const int nBins = 500;
            TH1::StatOverflows(true);
            TH1F hElectrons("hElectrons", "Number of electrons", nBins, 0., nBins);
            const int nEvents = 1000000;
            int echecs = 0;
            int division = int(nEvents/10);
            for (unsigned int i = 0; i < nEvents; ++i) {
                if (i % division == 0) cout << i << "/" << nEvents << "\n";
                // Initial coordinates of the photon.
                double x0 = width/2.;
                double y0 = depth/2.;
                double z0 = ddrift;
                double t0 = 0.;
                // Sample the photon energy, using the relative intensities according to XDB.
                const double r = 167. * RndmUniform();
                const double egamma = r < 100. ? 5898.8 : r < 150. ? 5887.6 : 6490.4;
                /* // translation
                 double egamma;
                 if ( r<100) egamma = 5898.8;
                 else if ( r<150 ) egamma = 5887.6;
                 else egamma = 6490.4;
                 */
                int ne = 0;
                track.TransportPhoton(x0, y0, z0, t0, egamma, 0., 0., -1, ne);
                if (ne < 3) echecs++;
            }
            hConversion->Fill(ddrift, 1-(double)echecs/nEvents);
            cout << "\n" << nEvents << " events simulated" << endl;
            cout << "\n" << echecs << " photons did not convert into electrons" << endl;
        }
        hConversion->Write();
        fOut->Close();
    }
    
    else {
        /*
        // Create a cylinder in which the x-rays can convert.
        // Diameter [cm]
        const double diameter = 7.8;
        // Half-length of the cylinder [cm].
        const double length = 10.;
        SolidTube tube(0.0, 0.0, 0.0, 0.0, 0.5 * diameter, length);
        
        // Combine gas and box to a simple geometry.
        geo->AddSolid(&tube, gas);
         */
        
        // Set up detector geometry
        GeometrySimple* geo = new GeometrySimple();
        // Set up the bulk
        SolidBox bulk = SolidBox(width/2., depth/2., ddrift/2., width/2., depth/2., ddrift/2.);
        geo->AddSolid(&bulk, gas);
        
        // Make a component
        ComponentConstant field;
        field.SetGeometry(geo);

        // Make a sensor.
        Sensor sensor;
        sensor.AddComponent(&field);

        // Use Heed for simulating the photon absorption.
        TrackHeed track;
        track.SetSensor(&sensor);
        //track.EnableElectricField();
        // Histogram
        const int nBins = 500;
        TH1::StatOverflows(true);
        TH1F hElectrons("hElectrons", "Number of electrons", nBins, 0., nBins);
        const int nEvents = 1000000;
        int division = int(nEvents/10);
        for (unsigned int i = 0; i < nEvents; ++i) {
            if (i % division == 0) cout << i << "/" << nEvents << "\n";
            // Initial coordinates of the photon.
            double x0 = width/2.;
            double y0 = depth/2.;
            double z0 = ddrift;
            /*
            double x0 = 0.;
            double y0 = 0.;
            double z0 = 0.;
             */
            double t0 = 0.;
            // Sample the photon energy, using the relative intensities according to XDB.
            const double r = 167. * RndmUniform();
            const double egamma = r < 100. ? 5898.8 : r < 150. ? 5887.6 : 6490.4;
            /* // translation
             double egamma;
             if ( r<100) egamma = 5898.8;
             else if ( r<150 ) egamma = 5887.6;
             else egamma = 6490.4;
             */
            int ne = 0;
            track.TransportPhoton(x0, y0, z0, t0, egamma, 0., 0., -1, ne);
            if (ne > 2) hElectrons.Fill(ne);
        }
    
        const bool drawSpectrum = false;
        if (drawSpectrum) {
            TCanvas c("c", "", 600, 600);
            c.cd();
            hElectrons.SetFillColor(kBlue + 2);
            hElectrons.SetLineColor(kBlue + 2);
            hElectrons.SetXTitle("# primary electrons");
            hElectrons.SetYTitle("# counts");
            hElectrons.Draw();
            c.SaveAs(Form("Figures/Fe55Spectrum-%s.pdf", gasName.c_str()));
            
            int nPrimary = hElectrons.GetMaximumBin();
            cout << "\nnPrimary = " << nPrimary << " in " << gasName << endl;
            
            cout << "\nCalculations yield:" << endl;
            cout << "nPrimary = 228 in Ar-iC4H10 (90/10)" << endl;
            cout << "nPrimary = 222 in Ar-CO2 (93/7)" << endl;
            cout << "nPrimary = 157 in Ne" << endl;
        }
    
        const bool saveResults = false;
        if (saveResults) {
        // Write the histograms to the TFile.
        const char* name = Form("rootFiles/%s/spectrum-Fe55.root", gasName.c_str());
        TFile* f = new TFile(name, "RECREATE");
        hElectrons.Write();
        f->Close();
        }
    }
    
    time_t t1 = time(NULL);
    PrintTime(t0, t1);
    
    //app.Run(true);
}


