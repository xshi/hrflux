// Author: Xin Shi <Xin.Shi@cern.ch>
// Based on draw.cc
// Created: [2013-05-10 Fri 14:57] 

#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TLegend.h>

using namespace std; 


TCanvas* myCanvas;

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
  TString figpath = "/afs/cern.ch/user/x/xshi/www/pxl/fig"; 
  TString figfile = Form("%s/effvsdc_%s.pdf", figpath.Data(), label.Data()); 
  TString rootfile = Form("%s/effvsdc_%s.root", figpath.Data(), label.Data()); 
  
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
  if ( label != "roc2_4") {
    cerr << "Please use: draw roc2_4" << endl; 
    return 1 ; 
  }

  draw(label); 

  gSystem->Exit(0);

  return 0 ;
}

#endif

