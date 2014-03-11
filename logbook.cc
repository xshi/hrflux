// Add contents in the logbook.root 
// Author: Xin Shi <Xin.Shi@cern.ch> 
// Created: <2013-03-07 Thu 13:34> 

#include <iostream>
#include "TTree.h"   
#include "TFile.h"   


using namespace std; 
// TGraph* test;

void logbook_v2(){

  TString infile =  "/home/pixel_dev/TB2012B_Data/logbook_v2.dat"; 
  TTree *T = new TTree("logbook", "logbook v2");
  Long64_t nlines = T->ReadFile(infile, "spill/I:x_mm/I:y_pos/I:PS/I:vthrcomp/I:HV/I:WBC/I:magslits_m/I:magslits_p/I:slits/I:xion/I");
  // Long64_t nlines = T->ReadFile(infile, "spill/I:x_mm/I:y_pos/I:PS/I:vthrcomp/I:HV/I:WBC/I:magslits_m/D:magslits_p/D:xion/I");
  
  printf(" found %lld points\n",nlines);
  
  TString outfile =  "/raid1/w/xshi/logbook/logbook_v2.root"; 

  TFile *f = new TFile(outfile, "recreate");
  
  T->Write();
  f->Close();

}
