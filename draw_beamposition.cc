// Author: Xin Shi <Xin.Shi@cern.ch>
// Based on draw.cc
// Created: [2013-05-28 Tue 10:50] 

#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TMath.h>
#include <TText.h>

// #include "fitspot.h"   


using namespace std; 


TCanvas* myCanvas;
TH1D *projX, *projY;
// TH1D *hposX, *hposY;
// TF2 *f2;
// TCanvas *myCanvas;
TF1 *myGausX, *myGausY; // , *myFunc2, *myFunc3; 

void set_root_style(int stat=1110, int grid=0){
  gROOT->Reset();

  gStyle->SetTitleFillColor(0) ; 
  gStyle->SetTitleBorderSize(0); 
    
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasDefX(0); 
  gStyle->SetCanvasDefY(0); 
  gStyle->SetFrameBorderMode(0); 
  gStyle->SetFrameBorderSize(1); 
  gStyle->SetFrameFillColor(0); 
  gStyle->SetFrameFillStyle(0); 
  gStyle->SetFrameLineColor(1); 
  gStyle->SetFrameLineStyle(1); 
  gStyle->SetFrameLineWidth(1); 

  gStyle->SetPadLeftMargin(0.10);  
  gStyle->SetPadRightMargin(0.05);  

  gStyle->SetLabelSize(0.03, "XYZ");  
  gStyle->SetTitleSize(0.04, "XYZ");  
  gStyle->SetTitleOffset(1.2, "Y");  

  gStyle->SetPadBorderMode(0);  
  gStyle->SetPadColor(0);  
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1); 
  gStyle->SetPadGridX(grid); 
  gStyle->SetPadGridY(grid); 

  gStyle->SetOptStat(stat); 
  gStyle->SetStatColor(0); 
  gStyle->SetStatBorderSize(1); 
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

double position(int runNumber, int nDetector, TString figfile, 
		double &posx, double &posy){

  // posx = 0, posy = 0; 
  // int runNumber = 16032;
  // int nDetector = 4; 
  // bool debug = false;
  // TString datadir="/home/pixel_dev/TB2012B_Data/data/"; 
  TString datadir="/Users/xshi/work/cms/pxl/dat/Run2011B/data/"; 


  // int bin = 3 ;
  // double posx, posy; 
  
  // TFile *myFile = TFile::Open(Form("%s%06d/%06d-clusters_bin%d.root", datadir.Data(), runNumber, runNumber, bin));

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
  TH2D *myHisto = (TH2D*)myDir->GetObjectChecked(Form("hitMap_d%d", nDetector), "TH2D");
  if (!myHisto) {
    std::cerr << "Could not find the histogram" << std::endl;
    return -1;
  }

  projX = myHisto->ProjectionX("projX");
  projY = myHisto->ProjectionY("projY");

  myGausX = new TF1("myGausX", "gaus;x [pixel unit]");
  // myGausY = new TF1("myGausY", "gaus");

  if (projX->GetEntries() == 0) return -1;
  if (projY->GetEntries() == 0) return -1;
    
  // if (debug) {
  // myCanvas = new TCanvas();
  //   myCanvas->cd();
  projX->Fit(myGausX);
  projX->SetTitle(Form("Run %d, ROC %d", runNumber, nDetector)); 
  myCanvas->Print(figfile); 

  // myCanvas = new TCanvas();
  // myCanvas->cd();
  // projY->Fit(myGausY);
  double xmin = 0;
  double xmax = 80; 
  int npar = 3; 

  TF1 *fy = new TF1("yfunc", AsyGausFixRatio, xmin, xmax, npar);

  fy->SetParNames("mean", "sigma1", "amplitude");
  fy->SetParameter("mean", 45);  
  fy->SetParameter("sigma1", 8);  
  fy->SetParameter("amplitude", 550);  
  projY->Fit("yfunc"); 
  projY->SetTitle(Form("Run %d, ROC %d", runNumber, nDetector)); 

  myCanvas->Print(figfile); 

  // } else {
  //   projX->Fit(myGausX, "QN");
  //   projY->Fit(myGausY, "QN");
  // }

  // double xMin, xMax;
  // xMin = projX->GetXaxis()->GetXmin();
  // xMax = projX->GetXaxis()->GetXmax();
  // double xFraction = myGausX->Integral(xMin, xMax)/myGausX->Integral(-1000, 1000);
  // xMin = projY->GetXaxis()->GetXmin();
  // xMax = projY->GetXaxis()->GetXmax();
  // double yFraction = myGausY->Integral(xMin, xMax)/myGausY->Integral(-1000, 1000);

  // double xyFraction = xFraction*yFraction; 

  posx = myGausX->GetParameter(1); 
  if (posx < 10) posx = projX->GetMean(); 
  posy = fy->GetParameter("mean"); 


  // printf("x = %4.2f, y = %4.2f\n", posx, posy);

  // if (posx < 0) posx = projX->GetMean();
  // if (posy < 0) posy = projY->GetMean(); 


  // return xyFraction;
  return 0; 
}



TEfficiency* drawEffvsDC(TTree* summary, TString figfile){
  const int NBINS = 26; 
  int xlow = 0; 
  int xup = 26; 

  int slits = 0;
  int slice = 0;
  int DC = 0;
  int nClusters = 0; 
  int nEvents = 0; 
  int nEmptyEvents = 0; 
  float eff_num = 0. ; 
  float eff_den = 0. ; 
  
  summary->SetBranchAddress("slits", &slits); 
  summary->SetBranchAddress("slice", &slice); 
  summary->SetBranchAddress("DC", &DC); 
  summary->SetBranchAddress("nClusters", &nClusters); 
  summary->SetBranchAddress("nEvents", &nEvents); 
  summary->SetBranchAddress("nEmptyEvents", &nEmptyEvents); 
  summary->SetBranchAddress("eff_num", &eff_num); 
  summary->SetBranchAddress("eff_den", &eff_den); 

  TH1D *h_pass = new TH1D("h_pass","pass", NBINS, xlow, xup);
  TH1D *h_total = new TH1D("h_total","total", NBINS, xlow, xup);

  Int_t nentries = (Int_t)summary->GetEntries();
  for (Int_t i=0;i<nentries;i++) {
    summary->GetEntry(i); 
    
    h_pass->Fill(DC, eff_num); 
    h_total->Fill(DC, eff_den); 
    
  }
  TEfficiency* pEff = new TEfficiency("eff", ";Double Column;Efficiency", NBINS, xlow, xup);

  for (Int_t bin=1; bin<=NBINS; bin++) {
    int iTotal = h_total->GetBinContent(bin); 
    int iPass = h_pass->GetBinContent(bin); 

    if (iPass  == 0 ) continue; 

    pEff->SetTotalEvents(bin, iTotal); 
    pEff->SetPassedEvents(bin, iPass); 
    
    // cout << "DC= " << bin << "pass = " << iPass << ", total = " << iTotal << 
    //   ", eff = " << float(iPass)/ iTotal << endl; 

  }

  delete h_pass; 
  delete h_total; 

  return pEff; 
}


void draw(TString label) {

  // TString figpath = "/afs/cern.ch/user/x/xshi/www/pxl/fig"; 

  TString figpath = "/Users/xshi/work/cms/pxl/fig"; 
  TString figfile = Form("%s/beamposition_%s.pdf", figpath.Data(), label.Data()); 
  TString rootfile = Form("%s/beamposition_%s.root", figpath.Data(), label.Data()); 

  cout << "label = " << label << endl; 

  double posx, posy; 
  
  int run = 16032; 
  
  myCanvas = new TCanvas("c", "c", 600, 600);
  set_root_style(); 
  myCanvas->UseCurrentStyle(); 
  myCanvas->Print(Form("%s[", figfile.Data()));


  
  TH2D *b2 = new TH2D("t2", ";x [pixel unit]; y [pixel unit]", 52, 0, 52, 80, 0, 80); 
  TText *a[8]; 

  for (int i=0; i < 8; i++) {
    position(run, i, figfile, posx, posy);
    // cout << "i = " << i << endl;  
    printf("x = %4.2f, y = %4.2f\n", posx, posy);
    b2->Fill(posx, posy); 
    b2->SetMarkerStyle(2);
    if (i == 1) {posx = posx+1;  posy = posy-2;}
    if (i == 2) {posx = posx+1;  posy = posy-.5;}
    if (i == 3) {posx = posx+1;  posy = posy-1;}
    if (i == 4) {posx = posx+1;  posy = posy-2;}
    if (i == 5) {posx = posx+1; }
    if (i == 7) {posx = posx-6;  posy = posy-2;}
    a[i] = new TText(posx, posy, Form("ROC %d", i)); 
    // a[i] = new TText(posx, posy, "ROC "); 


  }

  gStyle->SetOptStat(0); 
  myCanvas->UseCurrentStyle(); 
  b2->SetTitle(Form("Run %d", run)); 
  b2->SetMarkerStyle(kFullStar); 
  b2->SetMarkerSize(2); 
  b2->Draw(); 
  for (int i=0; i < 8; i++) {
    // a[i]->SetTextColor(i+1); 
    a[i]->SetTextSize(.03); 
    a[i]->Draw();     
  }


  myCanvas->Print(figfile); 
  
  myCanvas->Print(Form("%s]", figfile.Data()));

  // double  phi_D, phi_D_err, xposition_correction_factor; 
  // fitspot_v8(16032, 1, 4, phi_D, phi_D_err, xposition_correction_factor);

  // cout << phi_D << endl;
  /*
  double  posx = 0;
  double  posy = 0; 
  int runNumber = 16032;
  int nDetector = 4; 
  bool debug = false;
  TString datadir="/Users/xshi/work/cms/pxl/src/cc/"; 

  int bin = 1; 
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
  


  
  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double r = fitspot_v8(spill, i, nDetector, phi_D, phi_D_err, xposition_correction_factor);
  
  }

  TFile *_file0 = TFile::Open("/raid1/w/xshi/summary/summary_v4.root");
  TTree* summary0 = (TTree*) _file0->GetObjectChecked("summary", "TTree");

  TFile *_file1 = TFile::Open("/raid1/w/xshi/summary/summary_roc2.root");
  TTree* summary1 = (TTree*) _file1->GetObjectChecked("summary", "TTree");

  TFile *f = new TFile(rootfile, "RECREATE");

  TCanvas* myCanvas = new TCanvas("c", "c", 600, 600);
  set_root_style(); 
  myCanvas->UseCurrentStyle(); 
  myCanvas->Print(Form("%s[", figfile.Data()));

  TEfficiency* pEff = drawEffvsDC(summary0, figfile);
  TEfficiency* pEff1 = drawEffvsDC(summary1, figfile);

  TLegend* leg = new TLegend(0.1, 0.1, 0.2, .2);
  leg->AddEntry(pEff, "ROC4");
  leg->AddEntry(pEff1, "ROC2");
  pEff->SetLineColor(kBlue);
  pEff1->SetLineColor(kRed);
  pEff1->Draw("AP"); 
  pEff->Draw("same");
  leg->Draw(); 
  
  myCanvas->Print(figfile);

  myCanvas->Print(Form("%s]", figfile.Data()));
  f->Close();
  */ 
}



#ifndef __CINT__ 
#include <iostream>
#include <algorithm>

char* get_option(char ** begin, char ** end, const std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)  return *itr;
  return 0;
}

bool option_exists(char** begin, char** end, const std::string& option){
  return std::find(begin, end, option) != end;
}

void print_usage(){
  cerr << "Usage: draw version\n" << endl; 
}

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage() ;  
    return -1; 
  }

  TString label = argv[1]; 
  if ( label != "16032") {
    cerr << "Please use: draw_beamposition 16032" << endl; 
    return 1 ; 
  }

  draw(label); 

  gSystem->Exit(0);

  return 0 ;
}

#endif

