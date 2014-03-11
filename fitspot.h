#ifndef FITSPOT_H
#define FITSPOT_H

/* #include <TFile.h> */
/* #include <TCanvas.h> */
/* #include <TH1D.h> */
/* #include <TH2D.h> */
/* #include <TF2.h> */
/* #include <TAxis.h> */
/* #include <TDirectory.h> */
/* #include <TMath.h> */
/* #include <TSystem.h> */

// #include <TFitResultPtr.h>
/* #include <iostream> */
/* #include "TStyle.h" */
/* #include "TPaveText.h" */
/* #include <boost/format.hpp> */
/* #include <boost/filesystem/operations.hpp> */
/* #include <boost/filesystem.hpp> */

// using namespace std; 
// using namespace boost::filesystem;
// namespace fs = boost::filesystem;

// TFile *myFile;
/* TH2D *myHisto; */
/* TH1D *projX, *projY; */
/* TH1D *hposX, *hposY; */
/* TF2 *f2; */
/* TCanvas *myCanvas; */
/* TF1 *myGausX, *myGausY, *myFunc2, *myFunc3;  */

// TFitResultPtr rX; 
// TString datadir="/home/pixel_dev/TB2012B_Data/data/";  

static double _dummy_xyFraction;

double position(int bin, double &posx, double &posy);
double fitspot_v1(double &xyFraction=_dummy_xyFraction, 
		  int runNumber=16032, int nDetector=4, 
		  int verbose=1) ;

double fitspot_v2(double &xyFraction=_dummy_xyFraction, 
		  int runNumber=16032, int nDetector=4, 
		  int bin=-1, int verbose=1) ;

double fitspot_v3(double &xyFraction, double &xposition_correction_factor, 
		  int runNumber = 16032, int nDetector = 4, int bin=-1, int verbose=1); 
  
double fitspot_v4(double &xyFraction=_dummy_xyFraction, 
		  double &xposition_correction_factor=_dummy_xyFraction, 
		  int runNumber=16032, int nDetector=4, 
		  int bin=-1, int verbose=1) ;

double fitspot_v5(double &xyFraction=_dummy_xyFraction, 
		  double &xyFraction_err=_dummy_xyFraction, 
		  double &xposition_correction_factor=_dummy_xyFraction, 
		  int runNumber=16032, int nDetector=4, 
		  int bin=8, TString figname="", int verbose=3) ;

double fitspot_v8(int runNumber=15896, int bin=10, int nDetector=4, 
		  double &xyFraction=_dummy_xyFraction,
		  double &xyFraction_err=_dummy_xyFraction,
		  double &xposition_correction_factor=_dummy_xyFraction); 

void fitspot_bins_v8(int runNumber=15896, int nBins=24);

double peak_position_v5(int runNumber, int nDetector, int verbose=3); 
 
std::string check_and_join(std::string , std::string );

#endif

