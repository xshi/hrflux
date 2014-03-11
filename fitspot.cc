#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF2.h>
#include <TAxis.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TSystem.h>

// #include <TFitResultPtr.h>
#include <iostream>
#include "TStyle.h"
#include "TPaveText.h"
#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem.hpp>

#include "fitspot.h"


using namespace std; 
// using namespace boost::filesystem;
namespace fs = boost::filesystem;

// TFile *myFile;
TH2D *myHisto;
TH1D *projX, *projY;
TH1D *hposX, *hposY;
TF2 *f2;
TCanvas *myCanvas;
TF1 *myGausX, *myGausY, *myFunc2, *myFunc3; 

// TFitResultPtr rX; 
// TString datadir="/home/pixel_dev/TB2012B_Data/data/";  

// static double _dummy_xyFraction;

double macro(double &xyFraction, int runNumber = 16032, int nDetector = 4, 
	     int bin=-1,  bool debug = false) {

  TString filename; 
  TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 

  if (bin == -1)
    filename.Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber);

  else
    filename.Form("%s%06d/%06d-clusters_bin%d.root", datadir.Data(), runNumber, runNumber, bin); 

  // TFile *myFile = new TFile(); 
  // myFile->Open(filename, "READ"); 
  TFile * myFile = TFile::Open(filename);

  myFile = TFile::Open(Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber));
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }
  TDirectory* myDir = myFile->GetDirectory(Form("MyCMSPixelClusteringProcessor/detector_%d", nDetector));
  if (!myDir) {
    std::cerr << "Could not find the right directory" << std::endl;
    return -1;    
  }
  myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");
  myGausX = new TF1("myGausX", "gaus");
  myGausY = new TF1("myGausY", "gaus");
 
  if (projX->GetEntries() == 0) return -1;
  if (projY->GetEntries() == 0) return -1;

  if (debug) {
    myCanvas = new TCanvas();
    myCanvas->cd();
    // rX = projX->Fit(myGausX);
    projX->Fit(myGausX);
    myCanvas = new TCanvas();
    myCanvas->cd();
    projY->Fit(myGausY);
  } else {
    projX->Fit(myGausX, "QN");
    projY->Fit(myGausY, "QN");
  }

  double xMin, xMax;
  xMin = projX->GetXaxis()->GetXmin();
  xMax = projX->GetXaxis()->GetXmax();
  double xFraction = myGausX->Integral(xMin, xMax)/myGausX->Integral(-1000, 1000);
  xMin = projY->GetXaxis()->GetXmin();
  xMax = projY->GetXaxis()->GetXmax();
  double yFraction = myGausY->Integral(xMin, xMax)/myGausY->Integral(-1000, 1000);

  // double xyFraction = xFraction*yFraction; 
  xyFraction = xFraction*yFraction; 

  // myGausX->SetParameter(1, myGausX->GetParameter(1)*1.05); 
  // myGausY->SetParameter(1, myGausY->GetParameter(1)*1.05); 

  myGausX->SetParameter(1, projX->GetMean()); 
  myGausY->SetParameter(1, projY->GetMean()); 
  
  xFraction = myGausX->Integral(projX->GetXaxis()->GetXmin(), 
				projX->GetXaxis()->GetXmax())/myGausX->Integral(-1000, 1000);
  yFraction = myGausY->Integral(projY->GetXaxis()->GetXmin(), 
				projY->GetXaxis()->GetXmax())/myGausY->Integral(-1000, 1000);
  
  double xyFraction_new = xFraction*yFraction;

  double diff = (xyFraction_new - xyFraction)/xyFraction; 
  
  if (debug)
    cout << "diff = " << diff << endl; 

  myFile->Close(); 
  delete myFile; 
  return xyFraction;
}


void check_all(){
  int nBins = 24;
  // TFile *f = new TFile("check_all.root", "RECREATE");
  hposX = new TH1D("histo_pos","Position ;Time [s];  Position [mm]", nBins, 0, 12);
  hposY = new TH1D("histo_pos_Y","Y posion ;Time [s];Y position [pixel unit]", nBins, 0, 12);
  double posx, posy; 
  //double p = position(22, posx, posy);

  for (int i=1; i<=24; ++i) {
    // for (int i=22; i<=22; ++i) {
    double p = position(i, posx, posy);
    if (p == -1) {
      posx = 0;
      posy = 0;
    }
      
    // double p = i; 
    // hposX->Fill(p);
    // cout << i << ": posX = " << posx << ", posY = " << posy << endl; 
    hposX->SetBinContent(i, posx*0.15);
    hposY->SetBinContent(i, posy*0.1);
  }
  TCanvas *c1 = new TCanvas("c1", "hists with different scales",600,600);
  hposX->SetMaximum(80*0.1);
  hposY->SetMaximum(80*0.1);

  hposX->SetLineColor(kRed);
  hposX->Draw();
  c1->Update();

  hposY->SetLineColor(kBlue);
  hposY->Draw("same");

  // hposX->SaveAs("check_all.root");
  c1->SaveAs("posXY.pdf");
  // f->Close();

}


double position(int bin, double &posx, double &posy){

  posx = 0, posy = 0; 
  int runNumber = 16032;
  int nDetector = 4; 
  bool debug = false;
  TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 

  TFile *myFile = TFile::Open(Form("%s%06d/%06d-clusters_bin%d.root", datadir.Data(), runNumber, runNumber, bin));
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }
  TDirectory* myDir = myFile->GetDirectory(Form("MyCMSPixelClusteringProcessor/detector_%d", nDetector));
  if (!myDir) {
    std::cerr << "Could not find the right directory" << std::endl;
    return -1;    
  }
  myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");
  myGausX = new TF1("myGausX", "gaus");
  myGausY = new TF1("myGausY", "gaus");

  if (projX->GetEntries() == 0) return -1;
  if (projY->GetEntries() == 0) return -1;
    
  if (debug) {
    myCanvas = new TCanvas();
    myCanvas->cd();
    projX->Fit(myGausX);
    myCanvas = new TCanvas();
    myCanvas->cd();
    projY->Fit(myGausY);
  } else {
    projX->Fit(myGausX, "QN");
    projY->Fit(myGausY, "QN");
  }

  double xMin, xMax;
  xMin = projX->GetXaxis()->GetXmin();
  xMax = projX->GetXaxis()->GetXmax();
  double xFraction = myGausX->Integral(xMin, xMax)/myGausX->Integral(-1000, 1000);
  xMin = projY->GetXaxis()->GetXmin();
  xMax = projY->GetXaxis()->GetXmax();
  double yFraction = myGausY->Integral(xMin, xMax)/myGausY->Integral(-1000, 1000);

  double xyFraction = xFraction*yFraction; 

  printf("Bin %d : x= %4.2f, y=%4.2f, geo factor= %4.3f\n", bin,
	 myGausX->GetParameter(1), myGausY->GetParameter(1), xyFraction);
  
  posx =  myGausX->GetParameter(1); 
  posy =  myGausY->GetParameter(1); 

  if (posx < 0) posx = projX->GetMean();
  if (posy < 0) posy = projY->GetMean(); 

  return xyFraction;
}


double fitspot_v1(double &xyFraction, int runNumber, int nDetector, int verbose) {

  TString filename; 
  TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 

  // filename.Form("%s%06d/%06d-clusters_bin%d.root", datadir.Data(), runNumber, runNumber, bin); 
  filename.Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber); 

  TFile * myFile = TFile::Open(filename);
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }
  TDirectory* myDir = myFile->GetDirectory(Form("MyCMSPixelClusteringProcessor/detector_%d", nDetector));
  if (!myDir) {
    std::cerr << "Could not find the right directory" << std::endl;
    return -1;    
  }
  myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");
  myGausX = new TF1("myGausX", "gaus");
  myGausY = new TF1("myGausY", "gaus");
 
  if (projX->GetEntries() == 0) return -1;
  if (projY->GetEntries() == 0) return -1;

  if (verbose == 0){
    projX->Fit(myGausX, "QN");
    projY->Fit(myGausY, "QN");
  }
  else {
    projX->Fit(myGausX);
    projY->Fit(myGausY);
  }

  // if (verbose > 1) {
  //   myCanvas = new TCanvas();
  //   myCanvas->cd();
  //   rX = projX->Fit(myGausX);
  //   myCanvas = new TCanvas();
  //   myCanvas->cd();
  //   projY->Fit(myGausY);
  // } else {
  //   projX->Fit(myGausX, "QN");
  //   projY->Fit(myGausY, "QN");
  // }

  double xMin, xMax;
  xMin = projX->GetXaxis()->GetXmin();
  xMax = projX->GetXaxis()->GetXmax();
  double xFraction = myGausX->Integral(xMin, xMax)/myGausX->Integral(-1000, 1000);
  xMin = projY->GetXaxis()->GetXmin();
  xMax = projY->GetXaxis()->GetXmax();
  double yFraction = myGausY->Integral(xMin, xMax)/myGausY->Integral(-1000, 1000);

  // double xyFraction = xFraction*yFraction; 
  xyFraction = xFraction*yFraction; 

  // myGausX->SetParameter(1, myGausX->GetParameter(1)*1.05); 
  // myGausY->SetParameter(1, myGausY->GetParameter(1)*1.05); 

  myGausX->SetParameter(1, projX->GetMean()); 
  myGausY->SetParameter(1, projY->GetMean()); 
  
  xFraction = myGausX->Integral(projX->GetXaxis()->GetXmin(), 
				projX->GetXaxis()->GetXmax())/myGausX->Integral(-1000, 1000);
  yFraction = myGausY->Integral(projY->GetXaxis()->GetXmin(), 
				projY->GetXaxis()->GetXmax())/myGausY->Integral(-1000, 1000);
  
  double xyFraction_new = xFraction*yFraction;

  double diff = (xyFraction_new - xyFraction)/xyFraction; 
  
  if (verbose > 0)
    printf("xFraction = %f, yFraction = %f, xyFraction = %f, diff = %f \n", 
	   xFraction, yFraction, xyFraction, diff); 

  // myFile->Close(); 
  // delete myFile; 
  return xyFraction;
}


double fitspot_v2(double &xyFraction, int runNumber, int nDetector, 
	      int bin,  int verbose) {

  TString filename; 
  TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 

  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111111);

  if (bin == -1)
    filename.Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber);

  else
    filename.Form("%s%06d/%06d-clusters_bin%d.root", datadir.Data(), runNumber, runNumber, bin); 

  // TFile *myFile = new TFile(); 
  // myFile->Open(filename, "READ"); 
  TFile * myFile = TFile::Open(filename);
  // myFile = TFile::Open(Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber));
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }
  TDirectory* myDir = myFile->GetDirectory(Form("MyCMSPixelClusteringProcessor/detector_%d", nDetector));
  if (!myDir) {
    std::cerr << "Could not find the right directory" << std::endl;
    return -1;    
  }
  myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");
  myGausX = new TF1("myGausX", "gaus");
  myGausY = new TF1("myGausY", "gaus");
 
  if (projX->GetEntries() == 0) return -1;
  if (projY->GetEntries() == 0) return -1;

  if (verbose < 2){
    projX->Fit(myGausX, "QN");
    projY->Fit(myGausY, "QN");
  }
  else {
    projX->Fit(myGausX);
    projY->Fit(myGausY);
  }

  double xMin, xMax;
  xMin = projX->GetXaxis()->GetXmin();
  xMax = projX->GetXaxis()->GetXmax();
  double xFraction = myGausX->Integral(xMin, xMax)/myGausX->Integral(-1000, 1000);
  xMin = projY->GetXaxis()->GetXmin();
  xMax = projY->GetXaxis()->GetXmax();
  double yFraction = myGausY->Integral(xMin, xMax)/myGausY->Integral(-1000, 1000);

  // double xyFraction = xFraction*yFraction; 
  xyFraction = xFraction*yFraction; 

  // myGausX->SetParameter(1, myGausX->GetParameter(1)*1.05); 
  // myGausY->SetParameter(1, myGausY->GetParameter(1)*1.05); 

  myGausX->SetParameter(1, projX->GetMean()); 
  myGausY->SetParameter(1, projY->GetMean()); 
  
  xFraction = myGausX->Integral(projX->GetXaxis()->GetXmin(), 
				projX->GetXaxis()->GetXmax())/myGausX->Integral(-1000, 1000);
  yFraction = myGausY->Integral(projY->GetXaxis()->GetXmin(), 
				projY->GetXaxis()->GetXmax())/myGausY->Integral(-1000, 1000);
  
  double xyFraction_new = xFraction*yFraction;

  double xyFraction_err = (xyFraction_new - xyFraction)/xyFraction; 
  
  // if (debug)
  //   cout << "diff = " << diff << endl; 

  if (verbose > 0)
    printf("bin %d: xFrac = %f, yFrac = %f, xyFrac = %f, err = %f \n", 
	   bin, xFraction, yFraction, xyFraction, xyFraction_err); 

  // myFile->Close(); 
  // delete myFile; 
  // return xyFraction;
  return 0;
}

void fitspot_bins_v2(int runNumber=16032, int nBins=24) {
  int nDetector = 4; 
  for (int i=1; i<=nBins; ++i) {
    double t1; // , t2; 
    fitspot_v2(t1, runNumber, nDetector, i);
  }
}

double fitspot_v3(double &xyFraction, double &xposition_correction_factor, 
	      int runNumber, int nDetector, int bin, int verbose) {
  TString filename; 
  TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111111);

  // cout << "debug >>>> verbose = " << verbose << endl; 
  if (bin == -1)
    filename.Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber);
  else
    filename.Form("%s%06d/%06d-clusters_bin%d.root", datadir.Data(), runNumber, runNumber, bin); 

  TFile * myFile = TFile::Open(filename);
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }

  TDirectory* myDir = myFile->GetDirectory(Form("MyCMSPixelClusteringProcessor/detector_%d", nDetector));
  if (!myDir) {
    std::cerr << "Could not find the right directory" << std::endl;
    return -1;    
  }
  myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");

  // myFunc2 = new TF1("stefano2",
  // 		    "[0]*exp(2.59349*exp(-0.5*((-x+3.24273-[1])/7.07486)**2)+2.07765*exp(-0.5*((-x+9.33060e-01-[1])/2.24067)**2)-4.21847)",
  // 		    -40, 40);

  // myFunc3 = new TF1("stefano3",
  // 		    "[0]*exp(2.59349*exp(-0.5*((-x+3.24273-[1])/7.07486/0.150)**2)+2.07765*exp(-0.5*((-x+9.33060e-01-[1])/2.24067/0.150)**2)-4.21847)",
  // 		    -40, 40);

  myFunc3 = new TF1("stefano3", "[0]*exp(2.59349*exp(-0.5*((-x+3.24273/.15+[1]/.15)/7.07486*.15)**2)+2.07765*exp(-0.5*((-x+9.33060e-01/.15+[1]/.15)/2.24067*.15)**2)-4.21847)+[2]",-10, 20);
  //The default parameters to start the fit should be
  myFunc3->SetParameters(1, 0, 0);

  if (verbose == 0)
    projX->Fit(myFunc3, "QN");
  else
    projX->Fit(myFunc3);

  projX->GetXaxis()->SetTitle("X [pixel unit]");
  // projX->SetTitle(Form("Run %06d", runNumber));

  myFunc3->SetParameter(2,0);
  // myFunc3->SetNpx(100);
  double bx = myFunc3->Eval(24.6);
  // double bmax = myFunc3->Eval(myFunc3->GetMaximumX());
  double bmax = myFunc3->GetMaximum(0, 40);
  
  // cout << "Bin: " << bin << ", bx/bmax = " << bx/bmax
  //      << " , bx =" << bx << ", bmax = " << bmax << endl; 

  xposition_correction_factor = bx/bmax; 

  // myCanvas->Print(figname);
  double xFraction, yFraction; 

  if (myFunc3->Integral(-1000, 1000) == 0)
    xFraction = 0; 
  else
    xFraction = myFunc3->Integral(projX->GetXaxis()->GetXmin(), 
				  projX->GetXaxis()->GetXmax())/myFunc3->Integral(-1000, 1000);
  
  myGausY = new TF1("myGausY", "gaus");
  if (verbose == 0)
    projY->Fit(myGausY, "QN");
  else
    projY->Fit(myGausY);
  
  projY->GetXaxis()->SetTitle("Y [pixel unit]");
  // projY->SetTitle(Form("Run %06d", runNumber));

  if (myGausY->Integral(-1000, 1000) == 0)
    yFraction = 0;
  else
    yFraction = myGausY->Integral(projY->GetXaxis()->GetXmin(), 
				  projY->GetXaxis()->GetXmax())/myGausY->Integral(-1000, 1000);

  if (verbose > 1)
    printf("xFraction = %f, yFraction = %f ", xFraction, yFraction);

  xyFraction = xFraction*yFraction; 

  return 0; 
}

void fitspot_draw(TString figname="fitspot.pdf"){
  myCanvas = new TCanvas("fitspot", "fit beam spot", 600, 600);
  // fitspot_v3();  should be called seperately
  myCanvas->Print(Form("%s[", figname.Data()));
  projX->Draw(); 
  myCanvas->Print(figname);
  projY->Draw(); 
  myCanvas->Print(figname);
  myCanvas->Print(Form("%s]", figname.Data()));
}


void fitspot_bins_draw_v2(TString figname="fitspot_bins_v2.pdf", int runNumber=16032, int nBins=24) {
  myCanvas = new TCanvas("fitspot", "fit beam spot", 600, 600);
  int nDetector = 4; 
  
  myCanvas->Print(Form("%s[", figname.Data()));
  
  // nBins = 6; 
  
  int verbose = 2; // to draw the plots. 
  for (int i=1; i<=nBins; ++i) {
    // cout << bin << ": " << runNumber << endl; 
    double t1; // , t2; 
    double p = fitspot_v2(t1, runNumber, nDetector, i, verbose);
    if (p != 0) continue;
    projX->SetTitle(Form("Run %06d - Bin %d", runNumber, i));
    projX->Draw(); 
    myCanvas->Print(figname);
    
    projY->SetTitle(Form("Run %06d - Bin %d", runNumber, i));
    projY->Draw(); 
    myCanvas->Print(figname);
  }
  
  myCanvas->Print(Form("%s]", figname.Data()));
  
}

void fitspot_bins_draw_v3(TString figname="fitspot_bins_v3.pdf", int runNumber=16032, int nBins=24) {
  myCanvas = new TCanvas("fitspot", "fit beam spot", 600, 600);
  int nDetector = 4; 

  myCanvas->Print(Form("%s[", figname.Data()));
  
  // nBins = 4; 
  for (int i=1; i<=nBins; ++i) {
    // cout << bin << ": " << runNumber << endl; 
    double t1, t2; 
    double p = fitspot_v3(t1, t2, runNumber, nDetector, i);
    if (p != 0) continue;
    projX->SetTitle(Form("Run %06d - Bin %d", runNumber, i));
    projX->Draw(); 
    myCanvas->Print(figname);
    
    projY->SetTitle(Form("Run %06d - Bin %d", runNumber, i));
    projY->Draw(); 
    myCanvas->Print(figname);
  }
  
  myCanvas->Print(Form("%s]", figname.Data()));
  
}

double fitspot_v4(double &xyFraction, double &xposition_correction_factor, 
	      int runNumber, int nDetector, int bin, int verbose) {
  TString filename; 
  TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 
  // bin = 3; 
  // verbose = 3; 
  
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111111);

  // cout << "debug >>>> verbose = " << verbose << endl; 
  if (bin == -1)
    filename.Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber);
  else
    filename.Form("%s%06d/%06d-clusters_bin%d.root", datadir.Data(), runNumber, runNumber, bin); 

  TFile * myFile = TFile::Open(filename);
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }

  TDirectory* myDir = myFile->GetDirectory(Form("MyCMSPixelClusteringProcessor/detector_%d", nDetector));
  if (!myDir) {
    std::cerr << "Could not find the right directory" << std::endl;
    return -1;    
  }
  myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  // projX = myHisto->ProjectionX("projX", -40, 80);
  // projX->SetAxisRange(-40, 80); 
  
  projY = myHisto->ProjectionY("projY");

  // myFunc3 = new TF1("stefano3", "[0]*exp(2.59349*exp(-0.5*((-x+3.24273/.15+[1]/.15)/7.07486*.15)**2)+2.07765*exp(-0.5*((-x+9.33060e-01/.15+[1]/.15)/2.24067*.15)**2)-4.21847)+[2]",-10, 20);
  
  myFunc3 = new TF1("stefano3", "[0]*exp(2.59349*exp(-0.5*((-x+3.24273/.15+[1]/.15)/7.07486*.15)**2)+2.07765*exp(-0.5*((-x+9.33060e-01/.15+[1]/.15)/2.24067*.15)**2)-4.21847)+[2]",-40, 80);
  //The default parameters to start the fit should be
  myFunc3->SetParameters(1, 0, 0);

  if (verbose >= 3)
    projX->Fit(myFunc3);
  else
    projX->Fit(myFunc3, "QN");

  projX->GetXaxis()->SetTitle("X [pixel unit]");

  myFunc3->SetParameter(2,0);
  double bx = myFunc3->Eval(24.6);

  double bmax = myFunc3->GetMaximum(0, 40);
  
  xposition_correction_factor = bx/bmax; 
 
  double xFraction, yFraction; 

  if (myFunc3->Integral(-1000, 1000) == 0)
    xFraction = 0; 
  else
    // xFraction = myFunc3->Integral(projX->GetXaxis()->GetXmin(), 
    //				  projX->GetXaxis()->GetXmax())/myFunc3->Integral(-1000, 1000);
    xFraction = myFunc3->Integral(projX->GetXaxis()->GetXmin(), 
				  projX->GetXaxis()->GetXmax())/myFunc3->Integral(-20, 60);
  
  myGausY = new TF1("myGausY", "gaus");
  if (verbose >= 3)
    projY->Fit(myGausY);
  else
    projY->Fit(myGausY, "QN");
  
  projY->GetXaxis()->SetTitle("Y [pixel unit]");
  // projY->SetTitle(Form("Run %06d", runNumber));

  if (myGausY->Integral(-1000, 1000) == 0)
    yFraction = 0;
  else
    yFraction = myGausY->Integral(projY->GetXaxis()->GetXmin(), 
				  projY->GetXaxis()->GetXmax())/myGausY->Integral(-1000, 1000);

  if (verbose >=2)
    printf("xFraction = %f, yFraction = %f \n", xFraction, yFraction);

  xyFraction = xFraction*yFraction; 

  return 0; 
}


void fitspot_bins_v4(int runNumber=16032, int nBins=24) {
  int nDetector = 4; 
  int verbose = 2; 
  for (int i=1; i<=nBins; ++i) {
    double t1, t2; 
    fitspot_v4(t1, t2, runNumber, nDetector, i, verbose);
  }
}


void fitspot_bins_draw_v4(TString figname="~/www/pxl/fig/fitspot_bins_v4.pdf", 
			  int runNumber=16032, int nBins=24) {
  myCanvas = new TCanvas("fitspot", "fit beam spot", 600, 600);
  int nDetector = 4; 
  int verbose = 3; 
  myCanvas->Print(Form("%s[", figname.Data()));
  
  // nBins = 4; 
  for (int i=1; i<=nBins; ++i) {
    double t1, t2; 
    double p = fitspot_v4(t1, t2, runNumber, nDetector, i, verbose);
    if (p != 0) continue;
    projX->SetTitle(Form("Run %06d - Bin %d", runNumber, i));
    // projX->SetMinimum(0); 
    // projX->Draw(); 

    myFunc3->SetParameters(myFunc3->GetParameters());
    TH2F *h2tmp = new TH2F("h2tmp", "", 100, -40, 80, 100, -10., projX->GetMaximum()*1.2);
    h2tmp->Draw();
    myFunc3->SetLineColor(kBlue);
    myFunc3->Draw("same");
    projX->Draw("same");

    myCanvas->Print(figname);
    
    projY->SetTitle(Form("Run %06d - Bin %d", runNumber, i));
    projY->Draw(); 
    myCanvas->Print(figname);
  }
  
  myCanvas->Print(Form("%s]", figname.Data()));
  
}


double peak_position_v5(int runNumber, int nDetector, int verbose){
  double peak_position = 0; 
  // int runNumber = 16032;
  // int nDetector = 4; 
  TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 

  TFile *myFile = TFile::Open(Form("%s%06d/%06d-clusters.root", datadir.Data(), runNumber, runNumber));
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }

  TDirectory* myDir = myFile->GetDirectory(Form("MyCMSPixelClusteringProcessor/detector_%d", nDetector));
  if (!myDir) {
    std::cerr << "Could not find the right directory" << std::endl;
    return -1;    
  }
  myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  TH1D *projX_ = myHisto->ProjectionX("projX");
  // projX = myHisto->ProjectionX("projX", 10, 40);

  TF1 *myFunc = new TF1("stefano3", "[0]*exp(2.59349*exp(-0.5*((-x+3.24273/.15+[1]/.15)/7.07486*.15)**2)+2.07765*exp(-0.5*((-x+9.33060e-01/.15+[1]/.15)/2.24067*.15)**2)-4.21847)+[2]",-40, 80);
  myFunc->SetParameters(1, 0, 0);

  if (verbose == 0)
    projX_->Fit(myFunc, "QN");
  else
    projX_->Fit(myFunc);

  myFunc3->GetMaximumX(); 

  if ( TMath::Abs(myFunc->GetMaximumX() - projX->GetMaximumBin()) > 20){
    std::cerr << "Position fit is not correct for run: "  << runNumber << std::endl;
    delete projX_; 
    return -1;
  }

  double scintilator_B_off_set = 0.818096; // mm in x direction 
  // scintilator_B_off_set = scintilator_B_off_set / 0.15 ; // in pixel columns. 

  peak_position = myFunc->GetMaximumX() + scintilator_B_off_set / 0.15 ; 

  // projX->Clear(); 
  // myHisto->Clear();
  myFile->Close();
  delete myFunc; 
  delete projX_; 
  return peak_position; 
}

double fitspot_v5(double &xyFraction, double &xyFraction_err, double &xposition_correction_factor, 
		  int runNumber, int nDetector, int bin, TString figname,  int verbose) {
  // bin = 10; 
  // figname = "~/www/pxl/fig/test.pdf"; 
  xyFraction = 0;
  xposition_correction_factor = 1; 

  TString filename; 
  TString datadir="/raid1/w/ymtzeng/TIMER/HISTOGRAM_HitMapVSTime/";

  filename.Form("%s%06d_HitMapVSTime.root", datadir.Data(), runNumber);
    
  TFile * myFile = TFile::Open(filename);
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    return -1;
  }

  // myHisto = (TH2D*)myFile->GetObjectChecked(Form("ROC%dbin%d", nDetector, bin), "TH2D");
  gDirectory->GetObject(Form("ROC%dbin%d", nDetector, bin), myHisto); 

  if (!myHisto){
    std::cerr << "No object name found: " << Form("ROC%dbin%d", nDetector, bin) << std::endl;
    return -1;
  }

  if (myHisto->GetEntries() == 0) {
    printf("Empty entry for Run %d, bin %d. \n", runNumber, bin); 
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");

  projX->SetBinContent(1, projX->GetBinContent(1)/2.); 
  projX->SetBinContent(52, projX->GetBinContent(52)/2.); 
  projX->SetMinimum(0); 
  projY->SetMinimum(0); 

  //   myFunc3 = new TF1("stefano3", "[0]*exp(2.59349*exp(-0.5*((-x+3.24273/.15+[1]/.15)/7.07486*.15)**2)+2.07765*exp(-0.5*((-x+9.33060e-01/.15+[1]/.15)/2.24067*.15)**2)-4.21847)+[2]",-40, 80);
 
  myFunc3 = new TF1("stefano3", "[0]*exp(2.59349*exp(-0.5*((-x+3.24273/.15+[1]/.15)/7.07486*.15)**2)+2.07765*exp(-0.5*((-x+9.33060e-01/.15+[1]/.15)/2.24067*.15)**2)-4.21847)",-40, 80);
  //The default parameters to start the fit should be
  myFunc3->SetParameters(1, 0, 0);
  
  if (figname != "")
    projX->Fit(myFunc3);
  else
    projX->Fit(myFunc3, "QN");

  projX->GetXaxis()->SetTitle("X [pixel unit]");

  double xFraction=0; 
  double yFraction=0; 

  // myFunc3->SetParameter(2,0);
  // double peak_position = peak_position_v5(runNumber, nDetector, verbose-1);
  double peak_position = 24.6; 

  if (peak_position == -1)
    // Poisition fit is not good. Set the correctino factor 1. 
    return -1; 

  if (verbose > 2)
    printf("Beam spot peak position: %f [pixel unit].\n", peak_position); 
  
  // double bx = myFunc3->Eval(24.6);
  double bx = myFunc3->Eval(peak_position);
  
  //  double bmax = myFunc3->GetMaximum(0, 40);
  double bmax = myFunc3->GetMaximum(0, 52);
  xposition_correction_factor = bx/bmax; 
  
  
  xFraction = myFunc3->Integral(projX->GetXaxis()->GetXmin(), projX->GetXaxis()->GetXmax())
    /myFunc3->Integral(-300, 300);
  
  // double xFraction, xFraction_err, yFraction; 
  
  // if (myFunc3->Integral(-300, 300) == 0)
  //   xFraction = 0, xFraction_err = 0; 
  // else
  //   xFraction = myFunc3->Integral(projX->GetXaxis()->GetXmin(), projX->GetXaxis()->GetXmax())
  //     /myFunc3->Integral(-300, 300);
  
    // xFraction = myFunc3->IntegralAndError(projX->GetXaxis()->GetXmin(), 
    // 					  projX->GetXaxis()->GetXmax(), xFraction_err)
    //   /myFunc3->Integral(-300, 300);
  
  myGausY = new TF1("myGausY", "gaus");
  if (figname != "")
    projY->Fit(myGausY);
  else
    projY->Fit(myGausY, "QN");
  
  projY->GetXaxis()->SetTitle("Y [pixel unit]");
  // projY->SetTitle(Form("Run %06d", runNumber));

  if (myGausY->Integral(-1000, 1000) == 0)
    yFraction = 0;
  else
    yFraction = myGausY->Integral(projY->GetXaxis()->GetXmin(), 
				  projY->GetXaxis()->GetXmax())/myGausY->Integral(-1000, 1000);

  xyFraction = xFraction*yFraction; 

  // Get the xyFraction Error
  // myFunc3->SetParameter((myFunc3->GetParameter(0)+myFunc3->GetParError(0)),
  // 			(myFunc3->GetParameter(1)+myFunc3->GetParError(1))); 

  // myGausY->SetParameter(1, projY->GetMean()); 
  
  // double xFraction_new = myFunc3->Integral(projX->GetXaxis()->GetXmin(), projX->GetXaxis()->GetXmax())
  //   /myFunc3->Integral(-300, 300);
  // double yFraction_new = myGausY->Integral(projY->GetXaxis()->GetXmin(), projY->GetXaxis()->GetXmax())
  //   /myGausY->Integral(-1000, 1000);
  
  // // double xFrac_err = 
  
  // double xyFraction_new = xFraction_new*yFraction_new;

  // xyFraction_err = (xyFraction_new - xyFraction)/xyFraction; 
  xyFraction_err = 0.05; // use fixed value for this version. 
  
  // if (debug)
  //   cout << "diff = " << diff << endl; 

  if (verbose >= 2)
    printf("bin %d: xFrac = %f, yFrac = %f, xyFrac = %f, err = %f \n", 
	   bin, xFraction, yFraction, xyFraction, xyFraction_err); 
  


  if (figname != ""){
    TCanvas *c = new TCanvas("fitspot", "fit beam spot", 1200, 600);
    gStyle->SetOptStat(1111111);
    gStyle->SetOptFit(1111111);

    c->Divide(2); 
    c->cd(1);
    projX->SetTitle(Form("Run %06d - Bin %d", runNumber, bin));
    projX->Draw(); 
    c->cd(2);
    projY->SetTitle(Form("Run %06d - Bin %d", runNumber, bin));
    projY->Draw(); 
    c->Print(figname);
    c->Close(); 
  }    
  
  delete myFunc3;
  delete myGausY; 
  // myHisto->Clear(); // will crash if calling peak_position_v5 
  projX->Clear(); 
  myFile->Close();
  return 0; 
}


Double_t AsyGaus(Double_t *xx, Double_t *par){
  Double_t x = xx[0]; 
  
  Double_t mean = par[0];
  Double_t sigma1 = par[1]; 
  Double_t sigma2 = par[2];
  Double_t amplitude = par[3];

  if (x < mean) 
    return amplitude * TMath::Gaus(x, mean, sigma1); 
  else
    return amplitude * TMath::Gaus(x, mean, sigma2); 
}

Double_t AsyGausFixRatio(Double_t *xx, Double_t *par){
  Double_t x = xx[0]; 

  Double_t ratio = 12.8/23.5; 

  Double_t mean = par[0];
  Double_t sigma1 = par[1]; 
  Double_t amplitude = par[2];

  Double_t sigma2 = sigma1/ratio; 

  if (x < mean) 
    return amplitude * TMath::Gaus(x, mean, sigma1); 
  else
    return amplitude * TMath::Gaus(x, mean, sigma2); 
}

// double get_myHisto(int runNumber=15895, int bin=10, int nDetector=4){
TFile * get_myHisto(int runNumber=15895, int bin=10, int nDetector=4){
  TString filename; 
  TString datadir="/raid1/w/ymtzeng/TIMER/HISTOGRAM_HitMapVSTime/";

  filename.Form("%s%06d_HitMapVSTime.root", datadir.Data(), runNumber);
    
  TFile * myFile = TFile::Open(filename);
  if (!myFile) {
    std::cerr << "Cannot open file for run " << runNumber << std::endl;
    // return -1;
    return NULL;
  }

  gDirectory->GetObject(Form("ROC%dbin%d", nDetector, bin), myHisto); 

  if (!myHisto){
    std::cerr << "No object name found: " << Form("ROC%dbin%d", nDetector, bin) << std::endl;
    // return -1;
    return NULL;
  }

  if (myHisto->GetEntries() == 0) {
    printf("Empty entry for Run %d, bin %d. \n", runNumber, bin); 
    // return -1;
    return NULL;
  }
  // return 0; 
  return myFile; 
}

double get_xFraction(double &xposition_correction_factor=_dummy_xyFraction){
  projX = myHisto->ProjectionX("projX");
  projX->SetBinContent(1, projX->GetBinContent(1)/2.); 
  projX->SetBinContent(52, projX->GetBinContent(52)/2.); 
  projX->SetMinimum(0); 
 
  TF1 *fx = new TF1("xfunc", "[0]*exp(2.59349*exp(-0.5*((-x+3.24273/.15+[1]/.15)/7.07486*.15)**2)+2.07765*exp(-0.5*((-x+9.33060e-01/.15+[1]/.15)/2.24067*.15)**2)-4.21847)",-40, 80);
  fx->SetParameters(1, 0, 0);

  projX->Fit(fx);

  projX->GetXaxis()->SetTitle("X [pixel]");

  double xFraction=0; 

  if (fx->Integral(-300, 300) != 0)
    xFraction = fx->Integral(projX->GetXaxis()->GetXmin(), projX->GetXaxis()->GetXmax())
      /fx->Integral(-300, 300);

  double bx = fx->Eval(24.6);
  double bmax = fx->GetMaximum(0, 52);
  xposition_correction_factor = bx/bmax; 

  if (xFraction > 1) {
    cerr << "xFraction > 1 !!! :" << xFraction 
	 << "Set it to 0. " << endl; 
    xFraction = 0 ; 
  }

  return xFraction; 
}

double get_yFraction(){
  projY = myHisto->ProjectionY("projY");
  double xmin = 0;
  double xmax = 80; 
  int npar = 4; 

  TF1 *fy = new TF1("yfunc", AsyGaus, xmin, xmax, npar);

  fy->SetParNames("mean", "sigma1", "sigma2", "amplitude");
  fy->SetParameter("mean", 45);  
  fy->SetParameter("sigma1", 8);  
  fy->SetParameter("sigma2", 15);  
  fy->SetParameter("amplitude", 550);  
  projY->Fit("yfunc"); 
  projY->GetXaxis()->SetTitle("Y [pixel]");

  double yFraction=0; 
  if (fy->Integral(-300, 300) != 0)
    yFraction = fy->Integral(projY->GetXaxis()->GetXmin(), projY->GetXaxis()->GetXmax())
      /fy->Integral(-300, 300);

  return yFraction; 
}


double get_yFraction_v2(){
  projY = myHisto->ProjectionY("projY");
  double xmin = 0;
  double xmax = 80; 
  int npar = 3; 

  TF1 *fy = new TF1("yfunc", AsyGausFixRatio, xmin, xmax, npar);

  fy->SetParNames("mean", "sigma1", "amplitude");
  fy->SetParameter("mean", 45);  
  fy->SetParameter("sigma1", 8);  
  // fy->SetParameter("sigma2", 15);  
  fy->SetParameter("amplitude", 550);  
  projY->Fit("yfunc"); 
  projY->GetXaxis()->SetTitle("Y [pixel]");

  double yFraction=0; 

  if (fy->Integral(-300, 300) != 0)
    yFraction = fy->Integral(projY->GetXaxis()->GetXmin(), projY->GetXaxis()->GetXmax())
      /fy->Integral(-300, 300);
  
  // cout << "yFraction in get_yFraction: " << yFraction << endl; 

  if (yFraction > 1) {
    cerr << "yFraction > 1 !!! :" << yFraction 
	 << "Set it to 0. " << endl; 
    yFraction = 0 ; 
  }
  return yFraction; 
}


TPaveText * get_FitParameters(TF1 *f1, int npar){
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->SetBorderSize(1);
  pt->SetFillColor(0); 
  // pt->SetLabel("Fit Parameters");
  Char_t message[80];
  sprintf(message, "#chi^{2}/ndf = %.2f / %d", f1->GetChisquare(), f1->GetNDF());
  pt->AddText(message);  
  for (int i=0; i < npar; i++)
    {
      sprintf(message,"%s = %.2f #pm %.2f", f1->GetParName(i),
	      f1->GetParameter(i), f1->GetParError(i)); 
      pt->AddText(message);  
    }
  
  return pt;
}


double fitspot_v7(int runNumber, int bin=10, int nDetector=4){
  gStyle->SetOptFit(0);
  // gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(1);
  
  // int runNumber=15895;
  // int bin=10;
  // int nDetector=4; 

  Char_t message[80];
  //  get_myHisto(runNumber, bin, nDetector);   
  // if (get_myHisto(runNumber, bin, nDetector) == -1)
  if (!get_myHisto(runNumber, bin, nDetector))
    return -1; 

  double xFraction = get_xFraction(); 
  double yFraction = get_yFraction(); 
    
  TCanvas *c = new TCanvas("fitspot", "fit beam spot", 800, 800);
  c->UseCurrentStyle();
  c->Divide(2, 2); 
  c->cd(1);

  TF1 *fx= projX->GetFunction("xfunc");
  projX->SetTitle(Form("Run %06d - Bin %d", runNumber, bin));
  projX->Draw(); 

  c->cd(2); 
  TPaveText * pt = get_FitParameters(fx, 2); 
  sprintf(message, "x fraction = %.2f", xFraction); 
  pt->AddText(message); 
  pt->Draw();

  TF1 *fy= projY->GetFunction("yfunc");
  c->cd(3);
  projY->SetTitle(Form("Run %06d - Bin %d", runNumber, bin));
  projY->Draw(); 
  c->cd(4); 
  pt = get_FitParameters(fy, 4); 
  sprintf(message, "y fraction = %.2f", yFraction); 
  pt->AddText(message); 
  pt->Draw();

  TString figname;
  figname.Form("/raid1/w/xshi/flux/v7/%06d/fitspot_bin_%d.pdf", runNumber, bin); 
  c->Print(figname);  
  // delete c; 
  return 0; 
}

void fitspot_bins_v7(int runNumber=15895, int nBins=24) {
  int nDetector = 4; 
  // for (int i=0; i<=nBins; ++i) {
  for (int i=2; i<=nBins-2; ++i) {
    fitspot_v7(runNumber, i, nDetector);
  }
}

double fitspot_v8(int runNumber, int bin, int nDetector, 
		  double &xyFraction,
		  double &xyFraction_err,
		  double &xposition_correction_factor){
  int ver = 8; 
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(1);
 
  Char_t message[80];
  
  //if (get_myHisto(runNumber, bin, nDetector) == -1)
  // if (!get_myHisto(runNumber, bin, nDetector))
  //   return -1; 

  TFile * myFile = get_myHisto(runNumber, bin, nDetector); 
  if (!myFile)
    return -1; 

  double xFraction = get_xFraction(xposition_correction_factor); 
 
  double yFraction = get_yFraction_v2(); 

  TCanvas *c = new TCanvas("fitspot", "fit beam spot", 800, 800);
  c->UseCurrentStyle();
  c->Divide(2, 2); 
  c->cd(1);

  TF1 *fx= projX->GetFunction("xfunc");
  projX->SetTitle(Form("Run %06d - Bin %d", runNumber, bin));
  projX->Draw(); 

  c->cd(2); 
  TPaveText * pt = get_FitParameters(fx, 2); 
  sprintf(message, "x fraction = %.2f", xFraction); 
  pt->AddText(message); 
  pt->Draw();

  TF1 *fy= projY->GetFunction("yfunc");
  c->cd(3);
  projY->SetTitle(Form("Run %06d - Bin %d", runNumber, bin));
  projY->Draw(); 

  c->cd(4); 
  pt = get_FitParameters(fy, 3); 

  sprintf(message, "y fraction = %.2f", yFraction); 

  pt->AddText(message); 
  pt->Draw();

  string p = (boost::format("/raid1/w/xshi/flux/v%d/%06d") 
	      % ver % runNumber).str();
  string f = (boost::format("fitspot_bin_%d.pdf") % bin).str(); 

  TString outfile(check_and_join(p, f)); 
  c->Print(outfile);  
  delete c; 
  delete myFile; 

  xyFraction = xFraction*yFraction; 
  xyFraction_err = 0.05; // use fixed value for this version. 
  
  return 0; 
}


string check_and_join(string p, string f){
  if (!fs::exists(p)) {
    if (!fs::create_directory(p)) {
      cerr << "Couldn't create directory " << p << endl;
      return 0;
    }
  }
    
  fs::path path(p); 
  fs::path name(f); 
  fs::path file = path / name; 
  return file.string(); 
}


void fitspot_bins_v8(int runNumber, int nBins) {
  int nDetector = 4; 
  for (int i=0; i<=nBins; ++i) {
    // for (int i=2; i<=nBins-2; ++i) {
    fitspot_v8(runNumber, i, nDetector);
  }
}


// #ifndef __CINT__ 
#if !defined(__CINT__) && !defined(__LIBS__)


char* get_option(char ** begin, char ** end, const std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)  return *itr;
  return 0;
}

bool option_exists(char** begin, char** end, const std::string& option){
  return std::find(begin, end, option) != end;
}

void print_usage(){
  cerr << "Usage: fitspot version spill bin \n"
       << "Versions: \n" 
       << " v8 \t\tAsymmtric Gaussian with fixed ratio\n"
       << endl; 
}


int main(int argc, char** argv) {
  if (argc < 3) {
    print_usage() ;  
    return -1; 
  }

  TString version = argv[1]; 
  if (version == "v8" and argc == 4) {
      int spill = atoi(argv[2]); 
      int bin   = atoi(argv[3]); 
      fitspot_v8(spill, bin); 
  }

  if (version == "v8" and argc == 3) {
      int spill = atoi(argv[2]); 
      fitspot_bins_v8(spill); 
  }

  else{
    print_usage() ; 
    return -1;
  }

  gSystem->Exit(0);

  return 0 ;
}

#endif

