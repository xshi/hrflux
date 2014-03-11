// Author: Xin Shi <Xin.Shi@cern.ch>
// Based on draw_effvsdc_flux.cc 
// Created: [2013-05-16 Thu 14:32] 

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
#include <TProfile.h>
#include <math.h> 

using namespace std; 

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



TEfficiency* eff_vs_flux(TTree* summary, TString title, int set_slice){
  double AREA = 0.624/26; // cm^2     4160*150*100*10e-8
  double TIME = 25 * pow(10, -9) ; // 25ns 

  const int NBINS = 100; 
  int xlow = 0; 
  int xup = 1000; 

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
    
    if (slice != set_slice) continue; 

    int total_events = nEvents + nEmptyEvents; 
    double nhits = nClusters*1.0/total_events; 
    double phi = nhits/(AREA*TIME)/1000000. ; // MHz 
    
    h_pass->Fill(phi, eff_num); 
    h_total->Fill(phi, eff_den); 
    
  }
  TEfficiency* pEff = new TEfficiency("eff", Form("%s;Double Column;Efficiency", title.Data())
				      , NBINS, xlow, xup);

  for (Int_t bin=1; bin<=NBINS; bin++) {
    int iTotal = h_total->GetBinContent(bin); 
    int iPass = h_pass->GetBinContent(bin); 

    if (iPass  == 0 ) continue; 

    pEff->SetTotalEvents(bin, iTotal); 
    pEff->SetPassedEvents(bin, iPass); 
  }

  delete h_pass; 
  delete h_total; 

  return pEff; 
}



void draw_slice(int slice, TString title, TFile *_file0, TCanvas* myCanvas,
		TString label, TString figfile) {

   TTree* summary0; 
   title = Form("%s, slice %d", title.Data(), slice); 
  // eff vs flux 
  summary0 = (TTree*) _file0->GetObjectChecked("summary", "TTree");
  TEfficiency* pEff = eff_vs_flux(summary0, title, slice);
  pEff->SetLineColor(kBlue);
  summary0->Clear(); 
  
  TProfile* frame = new TProfile("frame",
				 Form("%s;Flux [MHz/cm^{2}];Efficiency",
				      title.Data()), 100, 0, 1000, 0., 1);
  // frame->GetYaxis()->SetRangeUser(0.4, 1.0);

  frame->Draw();  
  pEff->Draw("same");
 
  myCanvas->Print(figfile);
  myCanvas->Write(); 

  delete pEff; 
  delete frame; 
}




void draw(TString label) {
  TString figpath = "/afs/cern.ch/user/x/xshi/www/pxl/fig"; 
  TString figfile = Form("%s/effvsflux_slices_%s.pdf", figpath.Data(), label.Data()); 
  TString rootfile = Form("%s/effvsflux_slices_%s.root", figpath.Data(), label.Data()); 
  
  TString summaryfile; 
  TString title; 

  if (label == "roc4") {
    summaryfile = "/raid1/w/xshi/summary/summary_v4.root"; 
    title = "ROC 4 as DUT";  
    }
  if (label == "roc2") {
    summaryfile = "/raid1/w/xshi/summary/summary_roc2.root"; 
    title = "ROC 2 as DUT";  
  }

  TFile *_file0 = TFile::Open(summaryfile);

  TFile *f = new TFile(rootfile, "RECREATE");

  TCanvas* myCanvas = new TCanvas("c", "c", 600, 600);
  set_root_style(0); 
  myCanvas->UseCurrentStyle(); 
  myCanvas->Print(Form("%s[", figfile.Data()));

  for (int slice = 0; slice < 24; slice++) {
    // for (int slice = 4; slice < 5; slice++) {
    draw_slice(slice, title, _file0, myCanvas, label, figfile); 
  }
 
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
  if ( label != "roc4" and label != "roc2" ) {
    cerr << "Support labels: roc4, roc2" << endl; 
    return 1 ; 
  }

  draw(label); 

  gSystem->Exit(0);

  return 0 ;
}

#endif

