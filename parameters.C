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


// Set up detector geometry
const double pitch = 0.0025;    // cm
const double damp = 0.0128;
const double ddrift = 0.5;      // cm
const double dmylar = 3.;       // cm
const double radius = 0.0004;   // cm
const int periodicityNum = 5000;
double width = periodicityNum * pitch;
double depth = periodicityNum * pitch;

// gas parameters
const double rPenning = 0.51;    // of the order of 40%
const double lambdaPenning = 0.;

