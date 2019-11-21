#include <iostream>
#include <fstream>
#include <ctime>

#include <TApplication.h>
#include <TCanvas.h>

#include "parameters.C"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewCell.hh"

using namespace Garfield;

/* problems to be solved:
 - plot evrything with high Nperiod
 - save canvas (for now it's empty)
 */

int main(int argc, char * argv[]) {
    
    //______________________
    // variables
    const int modelNum = 4;
    //____________________
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    // Set up detector geometry
    GeometrySimple* geo = new GeometrySimple();
    
    //double meshHV = -550.;
    
        //Load geometry parameters
    double damp = 0., ddrift = 0., dmylar = 0., radius = 0., pitch = 0., width = 0., depth = 0.;
    int periodicityNum = 0;
    LoadParameters(modelNum, periodicityNum, damp, ddrift, dmylar, radius, pitch, width, depth);
    
    // Setup a vacuum medium
    MediumMagboltz* m = new MediumMagboltz();

    // Set up the bulk
    SolidBox bulk = SolidBox(width/2., depth/2., ddrift/2., width/2., depth/2., ddrift/2.);
    geo->AddSolid(&bulk, m);
    

    for (int i = 0; i < periodicityNum; ++i) {
        //std::cout << i << " creating tubes" << std::endl;
        SolidTube* tubeX = new SolidTube(width/2., pitch/2. + i*pitch, damp, radius, width/2, 1, 0, 0);
        SolidTube* tubeY = new SolidTube(pitch/2. + i*pitch, width/2., damp, radius, width/2, 0, 1, 0);
        geo->AddSolid(tubeX, m);
        geo->AddSolid(tubeY, m);
        //std::cout << i << " Tubes created" << std::endl;
    }
    

    // Create a viewer.
    ViewGeometry* view = new ViewGeometry();
    // Set the pointer to the geometry.
    view->SetGeometry(geo);
    
    bool drawGeom = true;
    if (drawGeom) {
        TCanvas* cd = new TCanvas();
        view->SetCanvas(cd);
        view->Plot();
        cd->SaveAs("Figures/DetectorGeometry.pdf");
    }
    
    
    app.Run(kTRUE);
    
}

