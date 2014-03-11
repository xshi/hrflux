// Scan the position of the scintilator 
// Author: Xin Shi <Xin.Shi@cern.ch> 
// Created: <2013-01-23 Wed 13:59> 

#include <iostream>
#include "TTree.h"   
#include "TFile.h"   
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TStyle.h>

TTree *raw, *tree;
TFile *fout;
TH1D *projX, *projY;
TH1D *hposX, *hposY;

TCanvas *myCanvas;
TF1* myGausX, * myGausY;
TH2D *myHisto;

double position(double &posx, double &posy, int runNumber, bool debug=false, 
		TCanvas *myCanvas=NULL, int xmm=0, int ymm=0);

void posscan() {
  raw = new TTree("raw", "Raw data");
  raw->ReadFile("posscan.dat");
  bool debug = true; 
  bool test = false; 
  /*
  raw->SetMarkerStyle(7);
  raw->SetMarkerColor(kBlue);
  raw->Draw("log(A/xion):x_mm", "y_pos>2700&&y_pos<2850&&x_mm!=-2");
  //raw->Draw("log(A/xion):x_mm", "x_mm<10&&x_mm>-10&&x_mm!=-2");
  raw->SetMarkerColor(kRed);
  raw->Draw("log(B/xion)+3.5:x_mm", "y_pos>2700&&y_pos<2850&&x_mm!=-2", "same");
  raw->SetMarkerColor(kGreen);
  raw->Draw("log(B/xion):x_mm", "y_pos>2700&&y_pos<2850&&x_mm!=-2", "same");
  */ 

  Int_t spill, x_mm, y_pos;
  Long64_t B; 
  raw->SetBranchAddress("spill", &spill);
  raw->SetBranchAddress("x_mm", &x_mm);
  raw->SetBranchAddress("y_pos", &y_pos);
  raw->SetBranchAddress("B", &B);

  double Cx, Cy; 
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("spill", &spill, "spill/I");
  tree->Branch("x_mm", &x_mm, "x_mm/I");
  tree->Branch("B", &B, "B/L"); 
  tree->Branch("Cx", &Cx, "Cx/D");

  if (debug){
    myCanvas = new TCanvas("posscan", "position scan", 640, 640);
    myCanvas->Print("position.pdf[");
  }

  for (int i=0; i<raw->GetEntries(); i++) {
    raw->GetEntry(i);
    if (y_pos > 2700 && y_pos < 2850 && x_mm != -2) {
    // if (y_pos > 2700 && y_pos < 2850 && x_mm == -3) {

      if (debug){
	position(Cx, Cy, spill, debug, myCanvas, x_mm, y_pos);
	// myCanvas->Print("position.pdf");
      }
      else
	position(Cx, Cy, spill);
	
      Cx = Cx*0.15; //convert to [mm]
      // Cx = 0.15; //convert to [mm]
      // Cy = 0.1; 
      if (debug) 
	printf("i= %d, Spill %d :x= %4.2f, y=%4.2f \n", i, spill, Cx, Cy);
      tree->Fill();
      if (test && i>50) break; 
    }
  }
  if (debug)
    myCanvas->Print("position.pdf]");
  
  if (test)
    fout = new TFile("posscan_test.root", "RECREATE");
  else
    fout = new TFile("posscan.root", "RECREATE");

  if (debug){
    tree->SetMarkerStyle(7);
    tree->Draw("Cx:x_mm");
  }

  raw->Write();
  tree->Write(); 
  fout->Close();
  // exit(0);
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
