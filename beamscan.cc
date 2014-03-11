// Scan the position of the beam spot
// Author: Xin Shi <Xin.Shi@cern.ch> 
// Created: <2013-01-24 Thu 23:34>  

#include <iostream>
#include "TTree.h"   
#include "TFile.h"   
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TStyle.h>
#include <TSystem.h>

TTree *raw, *tree, *t;
TH1D *projX, *projY;
TH1D *hbeamx, *hbeamy;
TFile *fout; 
TCanvas *myCanvas;
TF1* myGausX, * myGausY;
TH2D *myHisto;
TString dir = "/afs/cern.ch/user/x/xshi/work/cms/pxl/dat/Run2011B/"; 

double position(double &posx, double &posy, int runNumber, bool debug=false, 
		TCanvas *myCanvas=NULL, int xmm=0, int ymm=0);

void draw(TString figname="/afs/cern.ch/user/x/xshi/public/html/beamscan.pdf"); 

void write(); 

void beamscan() {
  if (!gSystem->AccessPathName(dir+"beamscan.root",kFileExists))
    write();
  draw();
}

void write() {
  raw = new TTree("raw", "Raw data");
  raw->ReadFile("/afs/cern.ch/user/x/xshi/work/cms/pxl/dat/Run2011B/runlist_v1.dat");
  // bool debug = true; 

  int run;
  double beamx, beamy; 
  
  raw->SetBranchAddress("run", &run);
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("run", &run, "run/I");
  tree->Branch("beamx", &beamx, "beamx/D");
  tree->Branch("beamy", &beamy, "beamy/D");
    
  for (int i=0; i<raw->GetEntries(); i++) {
    raw->GetEntry(i);
    // beamx = i*10; 
    position(beamx, beamy, run); 
    tree->Fill(); 
  }

  fout = new TFile("/afs/cern.ch/user/x/xshi/work/cms/pxl/dat/Run2011B/beamscan.root", "RECREATE");
  raw->Write();
  tree->Write(); 
  fout->Close();
  // delete fout; 
}


double position(double &posx, double &posy, int runNumber, bool debug, TCanvas *c,
		int xmm, int ymm){
  posx = 0, posy = 0; 
  int nDetector = 4; 
  TString datadir="/home/pixel_dev/TB2012B_Data/data/";
  // bool debug = true;

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

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");
  myGausX = new TF1("myGausX", "gaus");
  myGausY = new TF1("myGausY", "gaus");

  if (projX->GetEntries() == 0) return -1;
  if (projY->GetEntries() == 0) return -1;
    
  if (debug) {
    gStyle->SetOptStat(1111111);
    gStyle->SetOptFit(1111111);
    
    // myCanvas = new TCanvas("posscan", "position scan", 640, 640);
    // myCanvas->Print(Form("./fig/%06d-position.pdf[", runNumber));
    
    // myCanvas = new TCanvas();
    // myCanvas->cd();

    projX->Fit(myGausX);
    projX->GetXaxis()->SetTitle("X [pixel unit]");
    projX->SetTitle(Form("Run %06d - X pos = %d", runNumber, xmm));
    projX->Draw(); 

    // myCanvas->Update();
    // myCanvas->Print(Form("./fig/%06d-position.pdf", runNumber));
    // // myCanvas = new TCanvas();
    // myCanvas->cd();
    c->Print("position.pdf");

    projY->Fit(myGausY);
    projY->GetXaxis()->SetTitle("Y [pixel unit]");
    projY->SetTitle(Form("Run %06d - Y pos = %d", runNumber, ymm));
    projY->Draw(); 

    // myCanvas->Update();
    // myCanvas->Print(Form("./fig/%06d-position.pdf", runNumber));
    // myCanvas->Print(Form("./fig/%06d-position.pdf]", runNumber));
    c->Print("position.pdf");

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

  if (debug)
    printf("Run %d :x= %4.2f, y=%4.2f, geo factor= %4.3f\n", runNumber,
	   myGausX->GetParameter(1), myGausY->GetParameter(1), xyFraction);
  
  posx =  myGausX->GetParameter(1); 
  posy =  myGausY->GetParameter(1); 

  if (posx < 0) posx = projX->GetMean();
  if (posy < 0) posy = projY->GetMean(); 

  return xyFraction;
}


void draw(TString figname) {
  int run;
  double beamx, beamy; 

  TFile *f = new TFile("/afs/cern.ch/user/x/xshi/work/cms/pxl/dat/Run2011B/beamscan.root");
  t = (TTree*)f->Get("tree");
  
  t->SetBranchAddress("run", &run);
  t->SetBranchAddress("beamx", &beamx);
  t->SetBranchAddress("beamy", &beamy);
    
  int nBins = 100; 
  hbeamx = new TH1D("hbeamx","Beam spot X;X [mm]", nBins, 0, 5);

  for (int i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);
    hbeamx->Fill(beamx*0.15);  
  }
  
  myCanvas = new TCanvas("beamscan", "beam spot scan", 600, 600);
  gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(1111111);
  hbeamx->Draw();
  myCanvas->Print(figname);

}
