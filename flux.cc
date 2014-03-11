// Calculat the flux based on the trigger data
// Author: Xin Shi <Xin.Shi@cern.ch> 
// Created: <2013-01-15 Tue 17:24> 

#include <iostream>
#include "TTree.h"   
#include "TFile.h"   
#include "TH1.h"   
#include "TH1D.h"   
#include "TProfile.h"   
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include <TMath.h>
#include "TStyle.h"
#include <boost/format.hpp>


#include "fitspot.h"   


using namespace std; 
// TGraph* test;
TH1D *hphi, *hflux, *hphid; 
TProfile *hprof2, *hprof3;
TTree *tree, *T, *Tx; 
// TCanvas *c; 
// TFile *f; defined in fitspot.cc 

void flux_v5(Int_t spill=16032, Int_t XION=379392, Int_t nBins=24,
	     TString rootfile="", 
	     TString fitspotfig="", 
	     int verbose=3); 

 
void flux(Int_t spill=16032, Int_t XION=379392, Int_t nBins=24) { 
  // TString datadir="/home/pixel_dev/TB2012B_Data/data/";  defined in the fitspot.cc 
  int nDetector = 4;
  // double fraction = macro(spill, nDetector, datadir); 

  // cout << "fraction: " << fraction << endl; 
  // exit(0); 
  int XION_RATE = 6850;
  // int XION = 0;
  // double phi_D = macro(spill, nDetector, datadir); 
  // double phi_D_err =  phi_D*0.05; // Using a constant 0.05 for now.  
 
  // if (spill == 16032)
  //   {
  //     XION = 379392; 
  //     phi_D = 0.826; 
  //     phi_D_err = 0.826*0.05; 
  //   }
  // TString dir = gSystem->UnixPathName(datadir); 
  // TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  // dir.ReplaceAll("/data/","/trigger/");
  
  // f = new TFile(Form("trigger_spill%d.root", spill),"RECREATE");
  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat", 
				     spill),
				"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");
  printf(" found %lld points\n",nlines);

  T->Write();

  Int_t time;
  double phi; 
    
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("time", &time, "time/I");
  tree->Branch("phi", &phi, "phi/D");
  
  // test = new TGraph();

  int B_counts;
  int T0_counts;
  int C_counts;
  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);
  T->SetBranchAddress("C", &C_counts);

  int lastCounts;
  int lastClock;
  int B_counts_total; 

  int lastCounts_C;

  T->GetEntry(0);
  lastCounts=B_counts;
  lastClock=T0_counts;
  lastCounts_C=C_counts;

  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  

  // cout << "B counts total: " << B_counts_total << endl;

  double B_derivative;
  double C_derivative;

  hphi = new TH1D("hphi","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  hflux = new TH1D("hflux","Flux distribution;Time [s];Flux [MHz]", nBins, 0, 12);
  // hphi->Sumw2(); 
  // hphi = new TH1D("hphi2","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  hprof2 = new TProfile("hprof_freq","Rate versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hprof3 = new TProfile("hprof_freq_sent","Rate sent versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  double binWidth = hphi->GetXaxis()->GetBinWidth(1);
  // std::cout << "Bin width is exactly " << binWidth << std::endl;

  // int iPoints=0;
  int time_zero = -1; 

  for (Long64_t i=1; i<nlines; ++i) {
  // for (Long64_t i=300; i<350; ++i) {
    T->GetEntry(i);
    B_derivative = (B_counts-lastCounts)/double(T0_counts-lastClock);
    C_derivative = (C_counts-lastCounts_C)/double(T0_counts-lastClock);

    // test->SetPoint(iPoints++, T0_counts, B_derivative);

    time = T0_counts; 
    phi = B_derivative/double(B_counts_total);

    if (time_zero == -1 && C_counts > 0){
      time_zero = T0_counts; 
      // cout << i << "Time=" << time << ", B = " << B_counts << ", C = " <<  C_counts << endl; 
    }


    // hphi->Fill(time, (B_counts-lastCounts)/double(B_counts_total)/binWidth);
    // hphi->Fill(time, (B_counts-lastCounts));
    // hphi->Fill(time);

    if (time_zero == -1)
      hphi->Fill(-1, (B_counts-lastCounts));
    else 
      {
	// time = time-time_zero; 
	hphi->Fill(time-time_zero, (B_counts-lastCounts));
      }

    hprof2->Fill(time-time_zero, B_derivative);
    hprof3->Fill(time-time_zero, C_derivative);
    tree->Fill();
    // hflux->Fill(1); 

    lastCounts=B_counts;
    lastClock=T0_counts;
    lastCounts_C = C_counts; 
  }
  
  
  for (Long64_t i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double e = hphi->GetBinError(i) ;
    // double err = sqrt(x);

    // hphi->SetBinContent(i, x/binWidth/B_counts_total);
    // hphi->SetBinError(i, e/binWidth/B_counts_total);
    hphi->SetBinContent(i, x/B_counts_total);
    hphi->SetBinError(i, e/B_counts_total);
    // cout << i << ": " << x << " +/- " << e << 
    //   " err: " << err << "  diff: " << (e-err) << endl; 

    double flux, flux_err, phi_T, phi_T_err ;
    phi_T = hphi->GetBinContent(i); 
    phi_T_err = hphi->GetBinError(i);
 
    // double phi_D = macro(spill, nDetector, datadir); 
    //    double phi_D = macro(spill, nDetector, datadir, i); 
    double phi_D;
    // double p = macro(phi_D, spill, nDetector, i); 

    int verbose = 0; 
    double xposition_correction_factor; 
    double p = fitspot_v3(phi_D, xposition_correction_factor, spill, nDetector, i, verbose); 
    // double p = macro(phi_D, spill, nDetector, datadir, -1); 
    double phi_D_err =  phi_D*0.05; // Using a constant 0.05 for now.  
    
    cout << "phi_D = " << phi_D 
	 << ", xpos corr factor = " << xposition_correction_factor 
	 << endl; 

    if (p == -1)  {
      flux = 0;
      flux_err = 0; 
    }
    else {
      flux = phi_T/binWidth * XION * XION_RATE * phi_D; 
      flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
    }
    hflux->SetBinContent(i, flux);
    // else {
    //   flux = 0;
    //   flux_err = 0; 
    // }
    // hflux->Fill(0);
    // if (flux != 0) {
    //   hflux->SetBinContent(i, flux);
    // // cout << "defult error: " << hflux->GetBinError(i) << ", calc error: " <<  flux_err << endl; 
    //   hflux->SetBinError(i, flux_err);
    // } 
    // else
    //   hflux->Fill(0);
      
    // cout << i << "t0 + " << double(i)*0.5 << ": " << x << endl;  
    printf("t0 + %2.1f : %10.8f +/- %10.8f \n", double(i)*0.5, flux, flux_err);
  }
  // tree->SetDirectory(f);
 
  // TFile *fout = new TFile(Form("~/work/cms/pxl/dat/Run2012B/trigger/trigger_spill%d.root", spill),"RECREATE");
  
  printf("Total flux = %f", hflux->Integral()); 

  TFile *f = new TFile("test.root", "RECREATE");
  
  tree->Write();
  hphi->Write();
  hflux->Write();
  hprof2->Write();
  hprof3->Write();

  f->Close(); 
  // TCanvas *c = new TCanvas("canvas", "Saved plots", 600, 600);
  // tree->Draw("phi:time", "", "prof");
  // c->Write(); 
  // c->SaveAs(Form("trigger_spill%d.pdf", spill)); 

}

void flux_v1(Int_t spill=16032, Int_t nBins=24) {

  int XION_RATE = 6850;
  int XION = 0;
  double phi_D = 0;
  double phi_D_err = 0; 
 
  if (spill == 16032)
    {
      XION = 379392; 
      phi_D = 0.826; 
      phi_D_err = 0.826*0.05; 
    }
  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat",spill),
				"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");
  printf(" found %lld points\n",nlines);

  Int_t time;
  double phi; 
    
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("time", &time, "time/I");
  tree->Branch("phi", &phi, "phi/D");

  int B_counts;
  int T0_counts;
  int C_counts;
  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);
  T->SetBranchAddress("C", &C_counts);

  int lastCounts;
  int lastClock;
  int B_counts_total; 

  int lastCounts_C;

  T->GetEntry(0);
  lastCounts=B_counts;
  lastClock=T0_counts;
  lastCounts_C=C_counts;

  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  

  // cout << "B counts total: " << B_counts_total << endl;

  double B_derivative;
  double C_derivative;

  hphi = new TH1D("histo_phi","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  hflux = new TH1D("hflux","Flux distribution;Time [s];Flux [MHz]", nBins, 0, 12);
  // hphi->Sumw2(); 
  hprof2 = new TProfile("hprof_freq","Rate versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hprof3 = new TProfile("hprof_freq_sent","Rate sent versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  double binWidth = hphi->GetXaxis()->GetBinWidth(1);
  // std::cout << "Bin width is exactly " << binWidth << std::endl;

  int time_zero = -1; 

  for (Long64_t i=1; i<nlines; ++i) {
    T->GetEntry(i);
    B_derivative = (B_counts-lastCounts)/double(T0_counts-lastClock);
    C_derivative = (C_counts-lastCounts_C)/double(T0_counts-lastClock);

    time = T0_counts; 
    phi = B_derivative/double(B_counts_total);

    if (time_zero == -1 && C_counts > 0){
      time_zero = T0_counts; 
      // cout << i << "Time=" << time << ", B = " << B_counts << ", C = " <<  C_counts << endl; 
    }

    if (time_zero == -1)
      hphi->Fill(-1, (B_counts-lastCounts));
    else 
	hphi->Fill(time-time_zero, (B_counts-lastCounts));

    hprof2->Fill(time-time_zero, B_derivative);
    hprof3->Fill(time-time_zero, C_derivative);
    tree->Fill();

    lastCounts=B_counts;
    lastClock=T0_counts;
    lastCounts_C = C_counts; 
  }
  
  for (Long64_t i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double e = hphi->GetBinError(i) ;

    hphi->SetBinContent(i, x/B_counts_total);
    hphi->SetBinError(i, e/B_counts_total);

    double flux, flux_err, phi_T, phi_T_err ;
    phi_T = hphi->GetBinContent(i); 
    phi_T_err = hphi->GetBinError(i);
 
    flux = phi_T/binWidth * XION * XION_RATE * phi_D; 
    
      if (flux == 0)
	flux_err = 0;
      else
	flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
      //  flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
    
    printf("t0 + %2.1f : %10.8f +/- %10.8f \n", i*0.5, flux, flux_err);
    hflux->SetBinContent(i, flux);
    hflux->SetBinError(i, flux_err);
    
  }
  printf("Total flux = %.2f\n", hflux->Integral()); 

  TFile *f = new TFile("flux_v1.root", "RECREATE");
  
  tree->Write();
  hphi->Write();
  hflux->Write();
  hprof2->Write();
  hprof3->Write();

  f->Close(); 
}


void flux_draw_v1(TString figname="~/www/pxl/fig/flux_v1.pdf"){
  TFile *f = new TFile("~/www/pxl/fig/flux_v1.root");
  hflux = (TH1D*)f->Get("hflux");

  gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(1111111);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetStatX(0.4); // fX2NDC
  gStyle->SetStatY(0.8); // fY2NDC

  hflux -> GetYaxis() -> SetTitleOffset(1.3);
  hflux -> GetYaxis() -> SetLabelSize(0.04);
  hflux -> GetXaxis() -> SetLabelSize(0.04);
  hflux -> GetYaxis() -> SetTitleSize(0.04);
  hflux -> GetXaxis() -> SetTitleSize(0.04);
  
  TCanvas *c = new TCanvas("flux_v1", "flux distribution v1", 600, 600);
  c->UseCurrentStyle();
  hflux->Draw("E");
  c->Print(figname);
}


void flux_v2(Int_t spill=16032, Int_t XION=379392, Int_t nBins=24) { 
  int nDetector = 4;
  int XION_RATE = 6850;
  // int XION = 0;
  // double phi_D = macro(spill, nDetector, datadir); 
  // double phi_D_err =  phi_D*0.05; // Using a constant 0.05 for now.  
 
  // if (spill == 16032)
  //   {
  //     XION = 379392; 
  //     phi_D = 0.826; 
  //     phi_D_err = 0.826*0.05; 
  //   }

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat", 
				     spill),
				"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");
  printf(" found %lld points\n",nlines);

  //  T->Write();

  Int_t time;
  double phi; 
    
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("time", &time, "time/I");
  tree->Branch("phi", &phi, "phi/D");
  
  // test = new TGraph();

  int B_counts;
  int T0_counts;
  int C_counts;
  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);
  T->SetBranchAddress("C", &C_counts);

  int lastCounts;
  int lastClock;
  int B_counts_total; 

  int lastCounts_C;

  T->GetEntry(0);
  lastCounts=B_counts;
  lastClock=T0_counts;
  lastCounts_C=C_counts;

  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  

  // cout << "B counts total: " << B_counts_total << endl;

  double B_derivative;
  double C_derivative;

  hphi = new TH1D("hphi","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  hflux = new TH1D("hflux","Flux distribution;Time [s];Flux [MHz]", nBins, 0, 12);
  hprof2 = new TProfile("hprof_freq","Rate versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hprof3 = new TProfile("hprof_freq_sent","Rate sent versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hphid = new TH1D("hphid","Detector fraction distribution;Time [s];Fraction", nBins, 0, 12);
  double binWidth = hphi->GetXaxis()->GetBinWidth(1);
  // std::cout << "Bin width is exactly " << binWidth << std::endl;

  // int iPoints=0;
  int time_zero = -1; 

  for (Long64_t i=1; i<nlines; ++i) {
  // for (Long64_t i=300; i<350; ++i) {
    T->GetEntry(i);
    B_derivative = (B_counts-lastCounts)/double(T0_counts-lastClock);
    C_derivative = (C_counts-lastCounts_C)/double(T0_counts-lastClock);

    time = T0_counts; 
    phi = B_derivative/double(B_counts_total);

    if (time_zero == -1 && C_counts > 0){
      time_zero = T0_counts; 
      // cout << i << "Time=" << time << ", B = " << B_counts << ", C = " <<  C_counts << endl; 
    }

    if (time_zero == -1)
      hphi->Fill(-1, (B_counts-lastCounts));
    else 
      {
	// time = time-time_zero; 
	hphi->Fill(time-time_zero, (B_counts-lastCounts));
      }

    hprof2->Fill(time-time_zero, B_derivative);
    hprof3->Fill(time-time_zero, C_derivative);
    tree->Fill();
    // hflux->Fill(1); 

    lastCounts=B_counts;
    lastClock=T0_counts;
    lastCounts_C = C_counts; 
  }
  
  
  for (Long64_t i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double e = hphi->GetBinError(i) ;
    // double err = sqrt(x);

    // hphi->SetBinContent(i, x/binWidth/B_counts_total);
    // hphi->SetBinError(i, e/binWidth/B_counts_total);
    hphi->SetBinContent(i, x/B_counts_total);
    hphi->SetBinError(i, e/B_counts_total);
    // cout << i << ": " << x << " +/- " << e << 
    //   " err: " << err << "  diff: " << (e-err) << endl; 

    double flux, flux_err, phi_T, phi_T_err ;
    phi_T = hphi->GetBinContent(i); 
    phi_T_err = hphi->GetBinError(i);
 
    // double phi_D = macro(spill, nDetector, datadir); 
    //    double phi_D = macro(spill, nDetector, datadir, i); 
    double phi_D;
    // double p = macro(phi_D, spill, nDetector, i); 

    // int verbose = 0; 
    // double xposition_correction_factor; 
    double p = fitspot_v2(phi_D, spill, nDetector, i); 
    double phi_D_err =  phi_D*0.05; // Using a constant 0.05 for now.  
    
    hphid->SetBinContent(i, phi_D);
    hphid->SetBinError(i, phi_D_err);
    // cout << "phi_D = " << phi_D << endl; 

    if (p == -1)  {
      flux = 0;
      flux_err = 0; 
    }
    else {
      flux = phi_T/binWidth * XION * XION_RATE * phi_D; 
      if (flux == 0)
	flux_err = 0;
      else
	flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
    }
    hflux->SetBinContent(i, flux);
    hflux->SetBinError(i, flux_err);
    printf("t0 + %2.1f : %10.8f +/- %10.8f \n", double(i)*0.5, flux, flux_err);
  }

  printf("Total flux = %.2f\n", hflux->Integral()); 
  TFile *f = new TFile("flux_v2.root", "RECREATE");
  
  tree->Write();
  hphi->Write();
  hflux->Write();
  hprof2->Write();
  hprof3->Write();
  hphid->Write(); 
  f->Close(); 
}

void flux_draw_v2(TString figname="~/www/pxl/fig/flux_v2.pdf"){
  TFile *f = new TFile("~/www/pxl/fig/flux_v2.root");
  hflux = (TH1D*)f->Get("hflux");
  hphid = (TH1D*)f->Get("hphid");

  gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(1111111);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetStatX(0.4); // fX2NDC
  gStyle->SetStatY(0.8); // fY2NDC

  hflux -> GetYaxis() -> SetTitleOffset(1.3);
  
  TCanvas *c = new TCanvas("flux_v2", "flux distribution v2", 600, 600);
  c->UseCurrentStyle();

  c->Print(Form("%s[", figname.Data()));

  hflux->Draw("E");
  c->Print(figname);

  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.7);

  hphid->GetYaxis()->SetTitleOffset(1.3);
  hphid->Draw("E");
  c->Print(figname);

  c->Print(Form("%s]", figname.Data()));
}



void flux_v3(Int_t spill=16032, Int_t XION=379392, Int_t nBins=24, int verbose=1) { 
  int nDetector = 4;
  int XION_RATE = 6850;
  bool debug = true; 

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat", 
				     spill),
				"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");
  Int_t time;
  double phi; 
    
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("time", &time, "time/I");
  tree->Branch("phi", &phi, "phi/D");

  int B_counts;
  int T0_counts;
  int C_counts;
  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);
  T->SetBranchAddress("C", &C_counts);

  int lastCounts;
  int lastClock;
  int B_counts_total; 

  int lastCounts_C;

  T->GetEntry(0);
  lastCounts=B_counts;
  lastClock=T0_counts;
  lastCounts_C=C_counts;

  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  

  if (debug)
    printf("debug>>>> B counts total: %d\n", B_counts_total); 

  double B_derivative;
  double C_derivative;

  hphi = new TH1D("hphi","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  hflux = new TH1D("hflux","Flux distribution;Time [s];Flux [MHz]", nBins, 0, 12);

  hprof2 = new TProfile("hprof_freq","Rate versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hprof3 = new TProfile("hprof_freq_sent","Rate sent versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hphid = new TH1D("hphid","Detector fraction distribution;Time [s];Fraction", nBins, 0, 12);

  double binWidth = hphi->GetXaxis()->GetBinWidth(1);
  // std::cout << "Bin width is exactly " << binWidth << std::endl;

  // int iPoints=0;
  int time_zero = -1; 

  for (Long64_t i=1; i<nlines; ++i) {
    T->GetEntry(i);
    B_derivative = (B_counts-lastCounts)/double(T0_counts-lastClock);
    C_derivative = (C_counts-lastCounts_C)/double(T0_counts-lastClock);

    time = T0_counts; 
    phi = B_derivative/double(B_counts_total);

    if (time_zero == -1 && C_counts > 0){
      time_zero = T0_counts; 
      // cout << i << "Time=" << time << ", B = " << B_counts << ", C = " <<  C_counts << endl; 
    }

    if (time_zero == -1)
      hphi->Fill(-1, (B_counts-lastCounts));
    else 
      {
	// time = time-time_zero; 
	hphi->Fill(time-time_zero, (B_counts-lastCounts));
      }

    hprof2->Fill(time-time_zero, B_derivative);
    hprof3->Fill(time-time_zero, C_derivative);
    tree->Fill();

    lastCounts=B_counts;
    lastClock=T0_counts;
    lastCounts_C = C_counts; 
  }
  
  // Apply the xposition_correction_factor

  double B_counts_total_corrected = 0; 
  double xposition_correction_factor; 
  
  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    // double e = hphi->GetBinError(i) ;
    // B_counts_total_corrected += x;  // if placing here, will get the same totol as B counts. 
    double phi_D;

    fitspot_v3(phi_D, xposition_correction_factor, spill, nDetector, i, verbose-1); 

    // xposition_correction_factor = 1.0; // just test if this makes the total
				       // flux change. 

    hphi->SetBinContent(i, x/xposition_correction_factor);
    // if (debug)
    //   printf("debug >>> bin %d:  before correction: %f, after: %f \n",
    // 	     i, x, hphi->GetBinContent(i)); 
    B_counts_total_corrected += hphi->GetBinContent(i); // get the correct B_counts_total_corrected 
  }
  
  if (debug)
    printf("debug >>> B_counts_total_corrected = %f, B_counts_total = %d \n", B_counts_total_corrected, B_counts_total); 
  
  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double e = hphi->GetBinError(i) ;
    // hphi->SetBinContent(i, x/B_counts_total);
    // hphi->SetBinError(i, e/B_counts_total);
    hphi->SetBinContent(i, x/B_counts_total_corrected);
    hphi->SetBinError(i, e/B_counts_total_corrected);

    // if (debug)
    //   printf("debug >>> bin %d:  before correction: %f, after: %f \n\n",
    // 	     i, x/double(B_counts_total), hphi->GetBinContent(i)); 

    double flux, flux_err, phi_T, phi_T_err ;
    phi_T = hphi->GetBinContent(i); 
    phi_T_err = hphi->GetBinError(i);
 
    double phi_D;

    // verbose = 3; 
    // int verbose = 0; 
    // double xposition_correction_factor;  has been defined before. 
    double p = fitspot_v3(phi_D, xposition_correction_factor, spill, nDetector, i, verbose-1); 
    // double p = fitspot_v2(phi_D, spill, nDetector, i, verbose-1); 
    
    // printf("debug >>> phi_D = %f, xpos corr factor = %f\n",  phi_D, xposition_correction_factor);

    // phi_D = 0.826; 
    double phi_D_err =  phi_D*0.05; // Using a constant 0.05 for now.  
    hphid->SetBinContent(i, phi_D);
    hphid->SetBinError(i, phi_D_err);
    
    // if (debug)
    //   printf("debug >>> phi_D = %f, xpos corr factor = %f\n",  phi_D, xposition_correction_factor);

    if (p == -1)  {
      flux = 0;
      flux_err = 0; 
    }
    else {
      flux = phi_T/binWidth * XION * XION_RATE * phi_D; 
      if (flux == 0)
	flux_err = 0;
      else
	flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
    }
    hflux->SetBinContent(i, flux);
    hflux->SetBinError(i, flux_err);
    
    printf("t0 + %2.1f : %10.8f +/- %10.8f \n", double(i)*0.5, flux, flux_err);
  }

  printf("Total flux = %.2f\n", hflux->Integral()); 

  TFile *f = new TFile("~/www/pxl/fig/flux_v3.root", "RECREATE");
  T->Write(); 
  tree->Write();
  hphi->Write();
  hflux->Write();
  hprof2->Write();
  hprof3->Write();
  hphid->Write();

  f->Close(); 
  // delete tree;
  // delete hphi;
  // delete hflux;
  // delete hprof2;
  // delete hprof3; 
}

void flux_draw_v3(TString figname="~/www/pxl/fig/flux_v3.pdf"){
  TFile *f = new TFile("~/www/pxl/fig/flux_v3.root");
  hflux = (TH1D*)f->Get("hflux");
  hphid = (TH1D*)f->Get("hphid");

  gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(1111111);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetStatX(0.4); // fX2NDC
  gStyle->SetStatY(0.8); // fY2NDC

  hflux -> GetYaxis() -> SetTitleOffset(1.3);
  
  TCanvas *c = new TCanvas("flux_v3", "flux distribution v3", 600, 600);
  c->UseCurrentStyle();
  c->Print(Form("%s[", figname.Data()));


  hflux->Draw("E");
  c->Print(figname);

  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.7);

  hphid->GetYaxis()->SetTitleOffset(1.3);
  hphid->Draw("E");
  c->Print(figname);

  c->Print(Form("%s]", figname.Data()));
}


void flux_v4(Int_t spill=16032, Int_t XION=379392, Int_t nBins=24, int verbose=1) { 
  int nDetector = 4;
  int XION_RATE = 6850;
  bool debug = true; 

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat", 
				     spill),
				"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");
  Int_t time;
  double phi; 
    
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("time", &time, "time/I");
  tree->Branch("phi", &phi, "phi/D");

  int B_counts;
  int T0_counts;
  int C_counts;
  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);
  T->SetBranchAddress("C", &C_counts);

  int lastCounts;
  int lastClock;
  int B_counts_total; 

  int lastCounts_C;

  T->GetEntry(0);
  lastCounts=B_counts;
  lastClock=T0_counts;
  lastCounts_C=C_counts;

  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  

  if (debug)
    printf("debug>>>> B counts total: %d\n", B_counts_total); 

  double B_derivative;
  double C_derivative;

  hphi = new TH1D("hphi","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  hflux = new TH1D("hflux","Flux distribution;Time [s];Flux [MHz]", nBins, 0, 12);

  hprof2 = new TProfile("hprof_freq","Rate versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hprof3 = new TProfile("hprof_freq_sent","Rate sent versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hphid = new TH1D("hphid","Detector fraction distribution;Time [s];Fraction", nBins, 0, 12);

  double binWidth = hphi->GetXaxis()->GetBinWidth(1);
  // std::cout << "Bin width is exactly " << binWidth << std::endl;

  // int iPoints=0;
  int time_zero = -1; 

  for (Long64_t i=1; i<nlines; ++i) {
    T->GetEntry(i);
    B_derivative = (B_counts-lastCounts)/double(T0_counts-lastClock);
    C_derivative = (C_counts-lastCounts_C)/double(T0_counts-lastClock);

    time = T0_counts; 
    phi = B_derivative/double(B_counts_total);

    if (time_zero == -1 && C_counts > 0){
      time_zero = T0_counts; 
      // cout << i << "Time=" << time << ", B = " << B_counts << ", C = " <<  C_counts << endl; 
    }

    if (time_zero == -1)
      hphi->Fill(-1, (B_counts-lastCounts));
    else 
      {
	// time = time-time_zero; 
	hphi->Fill(time-time_zero, (B_counts-lastCounts));
      }

    hprof2->Fill(time-time_zero, B_derivative);
    hprof3->Fill(time-time_zero, C_derivative);
    tree->Fill();

    lastCounts=B_counts;
    lastClock=T0_counts;
    lastCounts_C = C_counts; 
  }
  
  // Apply the xposition_correction_factor

  double B_counts_total_corrected = 0; 
  double xposition_correction_factor; 
  double phi_D, phi_D_err;
  
  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    // double e = hphi->GetBinError(i) ;
    // B_counts_total_corrected += x;  // if placing here, will get the same totol as B counts. 

    fitspot_v4(phi_D, xposition_correction_factor, spill, nDetector, i, verbose-1); 

    // xposition_correction_factor = 1.0; // just test if this makes the total
				       // flux change. 
    phi_D_err =  phi_D*0.05; // Using a constant 0.05 for now.  
    hphid->SetBinContent(i, phi_D);
    hphid->SetBinError(i, phi_D_err);

    hphi->SetBinContent(i, x/xposition_correction_factor);
    // if (debug)
    //   printf("debug >>> bin %d:  before correction: %f, after: %f \n",
    // 	     i, x, hphi->GetBinContent(i)); 
    B_counts_total_corrected += hphi->GetBinContent(i); // get the correct B_counts_total_corrected 
  }
  
  if (debug)
    printf("debug >>> B_counts_total_corrected = %f, B_counts_total = %d \n", B_counts_total_corrected, B_counts_total); 
  
  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double e = hphi->GetBinError(i) ;
    // hphi->SetBinContent(i, x/B_counts_total);
    // hphi->SetBinError(i, e/B_counts_total);
    hphi->SetBinContent(i, x/B_counts_total_corrected);
    hphi->SetBinError(i, e/B_counts_total_corrected);

    // if (debug)
    //   printf("debug >>> bin %d:  before correction: %f, after: %f \n\n",
    // 	     i, x/double(B_counts_total), hphi->GetBinContent(i)); 

    double flux, flux_err, phi_T, phi_T_err ;
    phi_T = hphi->GetBinContent(i); 
    phi_T_err = hphi->GetBinError(i);
 
    // double phi_D;

    // verbose = 3; 
    // int verbose = 0; 
    // double xposition_correction_factor;  has been defined before. 
    //    double p = fitspot_v4(phi_D, xposition_correction_factor, spill, nDetector, i, verbose-1); 
    // double p = fitspot_v2(phi_D, spill, nDetector, i, verbose-1); 
    
    // printf("debug >>> phi_D = %f, xpos corr factor = %f\n",  phi_D, xposition_correction_factor);

    // phi_D = 0.826; 
    // double phi_D_err =  phi_D*0.05; // Using a constant 0.05 for now.  
    // hphid->SetBinContent(i, phi_D);
    // hphid->SetBinError(i, phi_D_err);
    
    phi_D = hphid->GetBinContent(i); 
    phi_D_err = hphid->GetBinError(i); 

    // if (debug)
    //   printf("debug >>> phi_D = %f, xpos corr factor = %f\n",  phi_D, xposition_correction_factor);

    double p = 0; 
    if (p == -1)  {
      flux = 0;
      flux_err = 0; 
    }
    else {
      flux = phi_T/binWidth * XION * XION_RATE * phi_D; 
      if (flux == 0)
	flux_err = 0;
      else
	flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
    }
    hflux->SetBinContent(i, flux);
    hflux->SetBinError(i, flux_err);
    
    printf("t0 + %2.1f : %10.8f +/- %10.8f \n", double(i)*0.5, flux, flux_err);
  }

  printf("Total flux = %.2f\n", hflux->Integral()); 

  TFile *f = new TFile("~/www/pxl/fig/flux_v4.root", "RECREATE");
  T->Write(); 
  tree->Write();
  hphi->Write();
  hflux->Write();
  hprof2->Write();
  hprof3->Write();
  hphid->Write();

  f->Close(); 
  // delete tree;
  // delete hphi;
  // delete hflux;
  // delete hprof2;
  // delete hprof3; 
}

void flux_draw_v4(TString figname="~/www/pxl/fig/flux_v4.pdf"){
  TFile *f = new TFile("~/www/pxl/fig/flux_v4.root");
  hflux = (TH1D*)f->Get("hflux");
  hphid = (TH1D*)f->Get("hphid");

  gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(1111111);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetStatX(0.4); // fX2NDC
  gStyle->SetStatY(0.8); // fY2NDC

  hflux -> GetYaxis() -> SetTitleOffset(1.3);
  
  TCanvas *c = new TCanvas("flux_v4", "flux distribution v4", 600, 600);
  c->UseCurrentStyle();
  c->Print(Form("%s[", figname.Data()));


  hflux->Draw("E");
  c->Print(figname);

  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.7);

  hphid->GetYaxis()->SetTitleOffset(1.3);
  hphid->Draw("E");
  c->Print(figname);

  c->Print(Form("%s]", figname.Data()));
}



void flux_v5(Int_t spill, Int_t XION, Int_t nBins, TString rootfile, 
	     TString fitspotfig, int verbose) { 
  int nDetector = 4;
  int XION_RATE = 6850;
  // verbose = 3; 

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat", 
				     spill),
				"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");
  Int_t time;
  double phi; 
  
  tree = new TTree("tree", "Analyzed tree");
  tree->Branch("time", &time, "time/I");
  tree->Branch("phi", &phi, "phi/D");

  int B_counts;
  int T0_counts;
  int C_counts;
  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);
  T->SetBranchAddress("C", &C_counts);

  int lastCounts;
  int lastClock;
  int B_counts_total; 

  int lastCounts_C;

  T->GetEntry(0);
  lastCounts=B_counts;
  lastClock=T0_counts;
  lastCounts_C=C_counts;

  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  

  if (verbose > 2)
    printf("B counts total: %d\n", B_counts_total); 

  double B_derivative;
  double C_derivative;

  hphi = new TH1D("hphi","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  hflux = new TH1D("hflux","Flux distribution;Time [s];Flux [MHz]", nBins, 0, 12);

  hprof2 = new TProfile("hprof_freq","Rate versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hprof3 = new TProfile("hprof_freq_sent","Rate sent versus time;Time [us];Rate [MHz]",nBins, 0, 12e6);
  hphid = new TH1D("hphid","Detector fraction distribution;Time [s];Fraction", nBins, 0, 12);

  double binWidth = hphi->GetXaxis()->GetBinWidth(1);
  int time_zero = -1; 

  for (Long64_t i=1; i<nlines; ++i) {
    T->GetEntry(i);
    B_derivative = (B_counts-lastCounts)/double(T0_counts-lastClock);
    C_derivative = (C_counts-lastCounts_C)/double(T0_counts-lastClock);

    time = T0_counts; 
    phi = B_derivative/double(B_counts_total);

    if (time_zero == -1 && C_counts > 0){
      time_zero = T0_counts; 
    }

    if (time_zero == -1)
      hphi->Fill(-1, (B_counts-lastCounts));
    else 
	hphi->Fill(time-time_zero, (B_counts-lastCounts));

    hprof2->Fill(time-time_zero, B_derivative);
    hprof3->Fill(time-time_zero, C_derivative);
    tree->Fill();

    lastCounts=B_counts;
    lastClock=T0_counts;
    lastCounts_C = C_counts; 
  }
  
  // Apply the xposition_correction_factor

  double B_counts_total_corrected = 0; 
  double xposition_correction_factor; 
  double phi_D, phi_D_err ;
  TString fitspotfig_bin; 

  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    if (fitspotfig != ""){
      fitspotfig_bin = Form("%s_bin_%d.pdf", fitspotfig.Data(), i);
      if (verbose > 2)
	printf("Saving fitspot as: %s ... \n", fitspotfig_bin.Data());  
      // continue;
    }
    
    double r = fitspot_v5(phi_D, phi_D_err, xposition_correction_factor, spill, nDetector, i,
			  fitspotfig_bin, verbose-1); 
    // Form("%s_bin_%d.pdf", fitspotfig_bin.Data(), i), 0); 
    if (r == -1){
      // problem in fitting the spot
      phi_D = 0;
      xposition_correction_factor = 1; 
    }
      
    phi_D_err = phi_D*0.05; // Using a constant 0.05 for now.  

    hphid->SetBinContent(i, phi_D);
    hphid->SetBinError(i, phi_D_err);

    hphi->SetBinContent(i, x/xposition_correction_factor);
    if (verbose > 2)
      printf("bin %d:  before correction: %f, after: %f \n",
     	     i, x, hphi->GetBinContent(i)); 
    B_counts_total_corrected += hphi->GetBinContent(i); // get the correct B_counts_total_corrected 
  }
  
  if (verbose > 2)
    printf("B_counts_total_corrected = %f, B_counts_total = %d \n", 
	   B_counts_total_corrected, B_counts_total); 
  

  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double e = hphi->GetBinError(i) ;
    hphi->SetBinContent(i, x/B_counts_total_corrected);
    hphi->SetBinError(i, e/B_counts_total_corrected);

    double flux, flux_err, phi_T, phi_T_err ;
    phi_T = hphi->GetBinContent(i); 
    phi_T_err = hphi->GetBinError(i);
 
    // verbose = 3; 
    // int verbose = 0; 
    // double xposition_correction_factor;  has been defined before. 
    // double p = fitspot_v5(phi_D, xposition_correction_factor, spill, nDetector, i, "", 0); 

    // double phi_D_err =  phi_D*0.05; 

    // hphid->SetBinContent(i, phi_D);
    // hphid->SetBinError(i, phi_D_err);

    phi_D = hphid->GetBinContent(i); 
    phi_D_err = hphid->GetBinError(i); 
    
    // if (debug)
    //   printf("debug >>> phi_D = %f, xpos corr factor = %f\n",  phi_D, xposition_correction_factor);
    
    // double p = 0; 
    // if (p == -1)  {
    //   flux = 0;
    //   flux_err = 0; 
    // }
    // else {
    flux = phi_T/binWidth * XION * XION_RATE * phi_D; 
    if (flux == 0)
      flux_err = 0;
    else
      flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
    //   }
    hflux->SetBinContent(i, flux);
    hflux->SetBinError(i, flux_err);
    
    if ( verbose > 0)
      printf("t0 + %2.1f : %.2f +/- %.2f \n", double(i)*0.5, flux, flux_err);
  }
  if ( verbose > 0)
    printf("Total flux = %.2f\n", hflux->Integral()); 

  if (rootfile != ""){
    // TString rootfile = Form("%s.root", fluxfig.Data()); 
    TFile *f = new TFile(rootfile, "RECREATE");
    T->Write(); 
    tree->Write();
    hphi->Write();
    hflux->Write();
    hprof2->Write();
    hprof3->Write();
    hphid->Write();

    f->Close(); 
    if (verbose > 0)
      printf("Saved to %s.\n", rootfile.Data()); 
    
    // delete f; No need if f->Close()
    delete T; 
    delete tree;
    delete hphi;
    delete hflux;
    delete hprof2;
    delete hprof3;
    delete hphid; 
  }
}

void flux_draw_v5(TString rootfile, TString figname=""){
  TFile *f = new TFile(rootfile);
 
  hflux = (TH1D*)f->Get("hflux");
  hphid = (TH1D*)f->Get("hphid");

  gStyle->SetOptFit(1111111);
  gStyle->SetOptStat(1111111);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetStatX(0.4); // fX2NDC
  gStyle->SetStatY(0.8); // fY2NDC

  hflux -> GetYaxis() -> SetTitleOffset(1.3);
  
  TCanvas *c = new TCanvas("flux", "flux distribution", 1200, 600);
  c->UseCurrentStyle();

  c->Divide(2);
  c->cd(1);
  hflux->Draw("E");

  c->cd(2);
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.7);
 
  hphid->GetYaxis()->SetTitleOffset(1.3);
  hphid->Draw("E");

  if (figname != ""){
    c->Print(figname);
    f->Close();
    delete c; 
  }

}


void flux_v6() { 
  TString rootfile = "/raid1/w/xshi/flux/flux_v6.root"; 
  
  int spills[] = {
    15881, 
    15882,
    15884,
    15885,
    15886,
    15887,
    15890,
    15891,
    15892,
    15893,
    15894,
    15895,
    15896,
    15897,
    15898,
    15899,
    15901,
    15905,
    15906,
    15907,
    15908,
    15911,
    15912,
    15913,
    15914 };
  
  int nspills = sizeof(spills) / sizeof(spills[0]); 
  Int_t nBins = 24; 

  TH1D *hphi_T = new TH1D("hphi_T","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  TH1D *hphi_T_B0 = new TH1D("hphi_T_B0","Phi(B0) versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  TH1D *hphi_T_B1 = new TH1D("hphi_T_B1","Phi(B1) versus time;Time [us];Trigger fraction", nBins, 0, 12e6);

  for (int i_spill = 0; i_spill < nspills; i_spill++){
    int spill = spills[i_spill]; 
    T = new TTree("raw", "raw data from ascii file");
    Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat", 
				       spill),
				  "A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");
    
    int B_counts, B0_counts, B1_counts;
    int T0_counts;
    int C_counts;
    T->SetBranchAddress("T0", &T0_counts); 
    T->SetBranchAddress("B", &B_counts);
    T->SetBranchAddress("C", &C_counts);
    T->SetBranchAddress("B0", &B0_counts);
    T->SetBranchAddress("B1", &B1_counts);
  
    int lastCounts_B, lastCounts_B0, lastCounts_B1;
    int B_counts_total, B0_counts_total, B1_counts_total; 

    T->GetEntry(0);
    lastCounts_B=B_counts;
    lastCounts_B0=B0_counts;
    lastCounts_B1=B1_counts;

    T-> GetEntry(nlines-1);
    B_counts_total = B_counts;  
    B0_counts_total = B0_counts;  
    B1_counts_total = B1_counts;  

    int time_zero = -1; 
   
    for (Long64_t i=1; i<nlines; ++i) {
      T->GetEntry(i);

      if (time_zero == -1 && C_counts > 0){
	time_zero = T0_counts; 
      }

      if (time_zero == -1){
	hphi_T->Fill(-1, (B_counts-lastCounts_B));
	hphi_T_B0->Fill(-1, (B0_counts-lastCounts_B0));
	hphi_T_B1->Fill(-1, (B1_counts-lastCounts_B1));
      } else { 
	hphi_T->Fill(T0_counts-time_zero, (B_counts-lastCounts_B));
	hphi_T_B0->Fill(T0_counts-time_zero, (B0_counts-lastCounts_B0));
	hphi_T_B1->Fill(T0_counts-time_zero, (B1_counts-lastCounts_B1));
      }
      
      lastCounts_B=B_counts;
      lastCounts_B0=B0_counts;
      lastCounts_B1=B1_counts;

    } // end loop for one fill

  } // end loop for all fills 

  // calculate the total trigger 
  double B_counts_total = 0; 
  double B0_counts_total = 0; 
  double B1_counts_total = 0; 
  for (int i=1; i<=nBins; ++i) {
    B_counts_total += hphi_T->GetBinContent(i); 
    B0_counts_total += hphi_T_B0->GetBinContent(i); 
    B1_counts_total += hphi_T_B1->GetBinContent(i); 
  }

  // set the bin content and error 
  
  double x, e;
  for (int i=1; i<=nBins; ++i) {
    x = hphi_T->GetBinContent(i) ;
    e = hphi_T->GetBinError(i) ;
    hphi_T->SetBinContent(i, x/B_counts_total);
    hphi_T->SetBinError(i, e/B_counts_total);

    x = hphi_T_B0->GetBinContent(i) ;
    e = hphi_T_B0->GetBinError(i) ;
    hphi_T_B0->SetBinContent(i, x/B0_counts_total);
    hphi_T_B0->SetBinError(i, e/B0_counts_total);

    x = hphi_T_B1->GetBinContent(i) ;
    e = hphi_T_B1->GetBinError(i) ;
    hphi_T_B1->SetBinContent(i, x/B1_counts_total);
    hphi_T_B1->SetBinError(i, e/B1_counts_total);
  }

  // hphi_T->Draw(); 

  TFile *f = new TFile(rootfile, "RECREATE");
  hphi_T->Write();
  hphi_T_B0->Write();
  hphi_T_B1->Write();
    
  f->Close(); 
  printf("Saved to %s.\n", rootfile.Data()); 
    
  
  delete hphi_T;

}

void fig_conic_a_vs_time(Int_t spill = 15896, Int_t nBins = 100){
  // Int_t spill = 15866; 
  // Int_t spill = 15896; 
  // Int_t nBins = 1000; 

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat",  spill),"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");

  int A_counts;
  int T0_counts;

  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("A", &A_counts);

  int lastCounts_A;
  int A_counts_total;

  TH1D *hphi_T = new TH1D("hphi_T",";Time [s];Coinc A", nBins, 0, 12);
  
  T->GetEntry(0);
  lastCounts_A=A_counts;
  
  T-> GetEntry(nlines-1);
  A_counts_total = A_counts;  
  
  for (Long64_t i=1; i<nlines; ++i) {
      T->GetEntry(i);
      hphi_T->Fill(T0_counts*1e-6, (A_counts-lastCounts_A));
      lastCounts_A=A_counts;
  }

  TCanvas *c = new TCanvas("flux", "flux distribution", 600, 600);
  // T->Draw("A:T0");
  hphi_T->Draw();
  c->Print("~/www/pxl/fig/conic_a_vs_time.pdf"); 
}


void fig_conic_a_vs_time_zoom(Int_t spill = 15896, Int_t nBins = 100){
  // Int_t spill = 15866; 
  // Int_t spill = 16032; 
  // Int_t nBins = 100; 

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat",  spill),"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");

  int A_counts;
  int T0_counts;

  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("A", &A_counts);

  int lastCounts_A;
  int A_counts_total;

  TH1D *hphi_T = new TH1D("hphi_T",";Time [s];Coinc A", nBins, 0, 2);
  
  T->GetEntry(0);
  lastCounts_A=A_counts;
  
  T-> GetEntry(nlines-1);
  A_counts_total = A_counts;  
  
  for (Long64_t i=1; i<nlines; ++i) {
      T->GetEntry(i);
      hphi_T->Fill(T0_counts*1e-6, (A_counts-lastCounts_A));
      lastCounts_A=A_counts;
  }

  TCanvas *c = new TCanvas("flux", "flux distribution", 600, 600);
  // T->Draw("A:T0");
  hphi_T->Draw();
  c->Print("~/www/pxl/fig/conic_a_vs_time_zoom.pdf"); 
}

void fig_conic_b_vs_time(Int_t spill = 15896, Int_t nBins = 100){
  // Int_t spill = 15866; 
  // Int_t spill = 15896; 
  // Int_t nBins = 1000; 

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat",  spill),"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");

  int B_counts;
  int T0_counts;

  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);

  int lastCounts_B;
  int B_counts_total;

  TH1D *hphi_T = new TH1D("hphi_T",";Time [s];Coinc B", nBins, 0, 12);
  
  T->GetEntry(0);
  lastCounts_B=B_counts;
  
  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  
  
  for (Long64_t i=1; i<nlines; ++i) {
      T->GetEntry(i);
      hphi_T->Fill(T0_counts*1e-6, (B_counts-lastCounts_B));
      lastCounts_B=B_counts;
  }

  TCanvas *c = new TCanvas("flux", "flux distribution", 600, 600);
  // T->Draw("A:T0");
  hphi_T->Draw();
  c->Print("~/www/pxl/fig/conic_b_vs_time.pdf"); 
}

void fig_conic_b_vs_time_zoom(Int_t spill = 15896, Int_t nBins = 100){
  // Int_t spill = 15866; 
  // Int_t spill = 15896; 
  // Int_t nBins = 1000; 

  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat",  spill),"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");

  int B_counts;
  int T0_counts;

  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);

  int lastCounts_B;
  int B_counts_total;

  TH1D *hphi_T = new TH1D("hphi_T",";Time [s];Coinc B", nBins, 0, 2);
  
  T->GetEntry(0);
  lastCounts_B=B_counts;
  
  T-> GetEntry(nlines-1);
  B_counts_total = B_counts;  
  
  for (Long64_t i=1; i<nlines; ++i) {
      T->GetEntry(i);
      hphi_T->Fill(T0_counts*1e-6, (B_counts-lastCounts_B));
      lastCounts_B=B_counts;
  }

  TCanvas *c = new TCanvas("flux", "flux distribution", 600, 600);
  // T->Draw("A:T0");
  hphi_T->Draw();
  c->Print("~/www/pxl/fig/conic_b_vs_time_zoom.pdf"); 
}


void get_phi_vs_time(Int_t spill=15895, Int_t nBins=24){
  T = new TTree("raw", "raw data from ascii file");
  Long64_t nlines = T->ReadFile(Form("/home/pixel_dev/TB2012B_Data/trigger/trigger_spill%d.dat", spill),
				"A0/I:A1/I:B0/I:B1/I:A/I:B/I:C/I:T0/I:T1/I");

  Int_t time;
  double phi; 

  int B_counts, T0_counts,  C_counts;
  T->SetBranchAddress("T0", &T0_counts); 
  T->SetBranchAddress("B", &B_counts);
  T->SetBranchAddress("C", &C_counts);
  
  int lastCounts, lastClock, lastCounts_C;
  int B_counts_total; 

  T->GetEntry(0);
  lastCounts=B_counts;
  lastClock=T0_counts;
  lastCounts_C=C_counts;

  T->GetEntry(nlines-1);
  B_counts_total = B_counts;  

  double B_derivative;

  hphi = new TH1D("hphi","Phi versus time;Time [us];Trigger fraction", nBins, 0, 12e6);
  int time_zero = -1; 

  for (Long64_t i=1; i<nlines; ++i) {
    T->GetEntry(i);
    B_derivative = (B_counts-lastCounts)/double(T0_counts-lastClock);

    time = T0_counts; 
    phi = B_derivative/double(B_counts_total);

    if (time_zero == -1 && C_counts > 0){
      time_zero = T0_counts; 
    }

    if (time_zero == -1)
      hphi->Fill(-1, (B_counts-lastCounts));
    else 
	hphi->Fill(time-time_zero, (B_counts-lastCounts));

    lastCounts=B_counts;
    lastClock=T0_counts;
    lastCounts_C = C_counts; 
  }
}


void apply_xposition_correction_factor(Int_t spill=15895, Int_t nBins=24, 
				       Int_t nDetector=4){

  double B_counts_total_corrected = 0; 
  double xposition_correction_factor; 
  double phi_D, phi_D_err ;
  TString fitspotfig_bin; 
  
  hphid = new TH1D("hphid","Detector fraction distribution;Time [s];Fraction", nBins, 0, 12);

  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double r = fitspot_v8(spill, i, nDetector, phi_D, phi_D_err, xposition_correction_factor);
  
   if (r == -1){
      phi_D = 0;
      xposition_correction_factor = 1; 
    }
      
    phi_D_err = phi_D*0.05; // Using a constant 0.05 for now.  

    hphid->SetBinContent(i, phi_D);
    hphid->SetBinError(i, phi_D_err);

    hphi->SetBinContent(i, x/xposition_correction_factor);
    if (true)
      printf("bin %d:  before correction: %f, after: %f \n",
     	     i, x, hphi->GetBinContent(i)); 
    B_counts_total_corrected += hphi->GetBinContent(i); 
  }

  for (int i=1; i<=nBins; ++i) {
    double x = hphi->GetBinContent(i) ;
    double e = hphi->GetBinError(i) ;
    hphi->SetBinContent(i, x/B_counts_total_corrected);
    hphi->SetBinError(i, e/B_counts_total_corrected);
    }
  
}


void calc_flux(Int_t nXion, Int_t nBins=24){
  const int XION_RATE = 6850;
  
  double flux, flux_err, phi_T, phi_T_err, phi_D, phi_D_err ;
  double binWidth = hphi->GetXaxis()->GetBinWidth(1);

  hflux = new TH1D("hflux","Flux distribution;Time [s];Flux [MHz]", nBins, 0, 12);
  
  for (int i=1; i<=nBins; ++i) {
    phi_T = hphi->GetBinContent(i); 
    phi_T_err = hphi->GetBinError(i);
    
    phi_D = hphid->GetBinContent(i); 
    phi_D_err = hphid->GetBinError(i); 
    
    flux = phi_T/binWidth * nXion * XION_RATE * phi_D; 
    if (flux == 0)
      flux_err = 0;
    else
      flux_err = sqrt(pow(flux, 2) *( pow(phi_T_err/phi_T, 2) + pow(phi_D_err/phi_D, 2))); 
   

    hflux->SetBinContent(i, flux);
    hflux->SetBinError(i, flux_err);
    
    printf("t0 + %2.1f : %.2f +/- %.2f \n", double(i)*0.5, flux, flux_err);
    
  }
 
  printf("Total flux = %.2f\n", hflux->Integral()); 
}

Int_t get_xion_number(Int_t spill){
  Int_t nXion = 0; 
  TString filename = "/raid1/w/xshi/logbook/logbook_v2.root"; 
  TFile* logbookFile = TFile::Open(filename);
  if (!logbookFile) 
    std::cerr << "Could not open logbook file" << std::endl;
 
  TTree* logbook; 
  gDirectory->GetObject("logbook", logbook); 
  if (!logbook) {
    std::cerr << "Could not find logbook tree" << std::endl;
    logbookFile->Close();
  }

  Int_t l_spill, l_xion;
  logbook->SetBranchAddress("spill",&l_spill);
  logbook->SetBranchAddress("xion",&l_xion);
 
  Long64_t nSpillEntries = logbook->GetEntries();
  for (Long64_t iSpill=0; iSpill<nSpillEntries; ++iSpill) { 
    logbook->GetEntry(iSpill);
    if (spill == l_spill ){ 
      // cout << l_spill << ": " << l_xion << endl; 
      nXion = l_xion; 
      break ; // continue; 
    }
  }
  
  logbookFile->Close();
  delete logbookFile; 
  return nXion; 
}


void save_to_root(TString rootfile){
  TFile *f = new TFile(rootfile, "RECREATE");
  T->Write(); 
  // tree->Write();
  hphi->Write();
  hflux->Write();
  // hprof2->Write();
  // hprof3->Write();
  hphid->Write();

  f->Close(); 
  printf("Saved to %s.\n", rootfile.Data()); 
    
  // delete f; No need if f->Close()
  delete T; 
  delete tree;
  delete hphi;
  delete hflux;
  delete hprof2;
  delete hprof3;
  delete hphid; 
}



void flux_spill(TString label, Int_t spill=15895){
  int ver = 8; 
  Int_t nDetector=4; 
  int nBins = 24; 
  if (label == "v8ROC2") {
    ver = 8; 
    nDetector = 2; 
  }


  get_phi_vs_time(spill); 
  apply_xposition_correction_factor(spill, nBins, nDetector); 
  
  Int_t nXion = get_xion_number(spill); 
  calc_flux(nXion); 
  
  string p = (boost::format("/raid1/w/xshi/flux/v%d/%06d") 
	      % ver % spill).str();
  string f = (boost::format("flux_%s.root") % label.Data()).str(); 

  TString outfile(check_and_join(p, f)); 
  
  save_to_root(outfile); 

}

vector<int> get_spills(int ver){
  TString filename; 
  if (ver == 8) 
    filename = "/raid1/w/xshi/logbook/logbook_v2.root"; 

  TFile* logbookFile = TFile::Open(filename);
  if (!logbookFile) 
    std::cerr << "Could not open logbook file" << std::endl;
 
  TTree* logbook; 
  gDirectory->GetObject("logbook", logbook); 
  if (!logbook) {
    std::cerr << "Could not find logbook tree" << std::endl;
    logbookFile->Close();
  }

  Int_t l_spill;
  vector<int> spills; 
  logbook->SetBranchAddress("spill",&l_spill);
 
  Long64_t nSpillEntries = logbook->GetEntries();
  for (Long64_t iSpill=0; iSpill<nSpillEntries; ++iSpill) { 
    logbook->GetEntry(iSpill);
    
    spills.push_back(l_spill); 
  }

  logbookFile->Close();
  delete logbookFile; 
  return spills; 
}
  

void flux(TString label){
  int ver = 8; 
  vector<int> spills = get_spills(ver); 

  for (vector<int>::iterator spill=spills.begin(); 
       spill != spills.end(); ++spill){
    // if ( *spill == 15880) continue; 

    cout << "processing " << *spill << " ..." << endl;
    
    flux_spill(label, *spill); 
    // break; 
  }

  cout << "Total spills: " << spills.size() << endl; 
}



#ifndef __CINT__ 
char* get_option(char ** begin, char ** end, const std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)  return *itr;
  return 0;
}

bool option_exists(char** begin, char** end, const std::string& option){
  return std::find(begin, end, option) != end;
}

void print_usage(){
  cerr << "Usage: flux label \n" << endl; 
}


int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage() ;  
    return -1; 
  }

  TString label = argv[1]; 
  if (label != "v8ROC2"){
    cerr << "please use ./flux v8ROC2! " << endl; 
    return -1; 
  }

  flux(label); 
  
    // // if (string(argv[1]) == "v8"){
    // if (argc == 2) 
    //   flux_v8(); 

    // if (argc == 3) {
    //   int spill = atoi(argv[2]); 
    //   // cout << "spill = " << spill << endl; 
    //   flux_v8_spill(spill); 
    // }
    
    // if (argc == 4) {
    //   int spill = atoi(argv[2]); 
    //   int nDetector = atoi(argv[3]); 
    //   // cout << "      // cout << "spill = " << spill << endl; 
    //   flux_v8_spill(spill); 
    // }
    
    

  // }
  
  // else if (string(argv[1]) == "fitspot"){
  //   TString version = argv[2]; 
  //   if (version == "v8" and argc == 5) {
  //     int spill = atoi(argv[3]); 
  //     int bin   = atoi(argv[4]); 
  //     fitspot_v8(spill, bin); 
  //   }
    
  //   if (version == "v8" and argc == 4) {
  //     int spill = atoi(argv[3]); 
  //     fitspot_bins_v8(spill); 
  //   }
    
  //   else{
  //     print_usage() ; 
  //     return -1;
  //   }
    
  // }
    
  // else{
  //   print_usage() ; 
  //   return -1;
  // }

  gSystem->Exit(0);

  return 0 ;
}

#endif

