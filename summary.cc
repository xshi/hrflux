// Calculat the flux based on the trigger data
// Author: Xin Shi <Xin.Shi@cern.ch>
// Based on Stefano Mersi's buildTree.cpp 
// Created: <2013-02-20 Wed 09:50> 


#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TTree.h>

#include <vector>
#include <iostream>

using namespace std; 

int minValue(int a, int b) { return (a<b)?a:b; }

TTree* get_logbook(TString label) {
  TTree* logbook = NULL; 
  TFile* logbookFile = NULL;

  TString filename;
  TString datadir="/raid1/w/xshi/logbook"; 
  filename.Form("%s/logbook_%s.root", datadir.Data(), label.Data()); 
  logbookFile = TFile::Open(filename);
  if (!logbookFile) {
    std::cerr << "Could not open logbook file" << std::endl;
    return NULL; 
  }

  gDirectory->GetObject("logbook", logbook); 

  if (!logbook){
    std::cerr << "No object name logbook! " << std::endl;
    return NULL;
  }
  return logbook; 
}


void summary(TString label) {
  TString outfile = Form("/raid1/w/xshi/summary/summary_%s.root", label.Data()); 
  
  TTree* logbook = get_logbook("v2"); 

  if (!logbook) return; 

  TTree* summary = new TTree("summary", "summary");
  summary->SetDirectory(0); 

  // Declaration of leaves types
  Int_t    spill;
  Int_t    slice;
  Int_t    DC; // Double Column 
  Int_t    x_mm;
  Int_t    y_pos;
  Int_t    PS;
  Int_t    vthrcomp;
  Int_t    HV;
  Int_t    WBC;
  Int_t    slits;
  Int_t    nHitPixels;
  Int_t    nClusters;
  Int_t    nEvents;
  Int_t    nEmptyEvents;
  Float_t  eff_num;
  Float_t  eff_den;

 // Set branch addresses.
  summary->Branch("spill", &spill, "spill/I");
  summary->Branch("slice", &slice, "slice/I");
  summary->Branch("DC", &DC, "DC/I");
  summary->Branch("x_mm", &x_mm, "x_mm/I");
  summary->Branch("y_pos", &y_pos, "y_pos/I");
  summary->Branch("PS", &PS, "PS/I");
  summary->Branch("vthrcomp", &vthrcomp, "vthrcomp/I");
  summary->Branch("HV", &HV, "HV/I");
  summary->Branch("WBC", &WBC, "WBC/I");
  summary->Branch("slits", &slits, "slits/I");
  summary->Branch("nHitPixels", &nHitPixels, "nHitPixels/I");
  summary->Branch("nClusters", &nClusters, "nClusters/I");
  summary->Branch("nEvents", &nEvents, "nEvents/I");
  summary->Branch("nEmptyEvents", &nEmptyEvents, "nEmptyEvents/I");
  summary->Branch("eff_num", &eff_num, "eff_num/F");
  summary->Branch("eff_den", &eff_den, "eff_den/F");

  // Set branch addresses for the input logbook TTree
  logbook->SetBranchAddress("spill",&spill);
  logbook->SetBranchAddress("x_mm",&x_mm);
  logbook->SetBranchAddress("y_pos",&y_pos);
  logbook->SetBranchAddress("PS",&PS);
  logbook->SetBranchAddress("vthrcomp",&vthrcomp);
  logbook->SetBranchAddress("HV",&HV);
  logbook->SetBranchAddress("WBC",&WBC);
  logbook->SetBranchAddress("slits",&slits);

  string thisDir = gSystem->pwd(); 
  Long64_t nSpillEntries = logbook->GetEntries();

  TString effdir = "/raid1/w/ymtzeng/rslu/efficiencyROC2"; 
  TString hitdir = "/raid1/w/ymtzeng/TIMER/HISTOGRAM_HitMapVSTime";

  for (Long64_t iSpill=0; iSpill<nSpillEntries; ++iSpill) { 
    
    logbook->GetEntry(iSpill);
    TFile* pixFile = TFile::Open(Form("%s/%06d/NumberOfPixelHitVStime.root", 
				      effdir.Data(), spill));
    if (! pixFile) 
      cout << Form("%s/%06d/NumberOfPixelHitVStime.root", 
		   effdir.Data(), spill) << endl; 


    TFile* cluFile; 
 
    cluFile = TFile::Open(Form("%s/%06d/output.root", effdir.Data(), spill));
    TFile* evtFile = TFile::Open(Form("%s/%06d/NumberOfEventVStime.root", 
				      effdir.Data(), spill));

    TFile* effFile = TFile::Open(Form("%s/%06d/effVStime.root", 
				      effdir.Data(), spill));
      
    TFile* hitFile = TFile::Open(Form("%s/%06d_HitMapVSTime.root",
				      hitdir.Data(), spill));
    
    if (!(pixFile&&cluFile&&evtFile&&effFile&&hitFile)) {
      std::cout << Form("One of the required files for spill %d is missing", spill)
		<< std::endl;
      continue;
    }
    
    TH1F* pixHisto = (TH1F*) pixFile->GetObjectChecked("h_numh", "TH1F");
    TH1F* evtHisto; 
    
    TH1D* cluHisto1;
    TH1D* nDoubleColumnsHisto; 
    TH2F * colvscluHisto; 

    TH2F * effNumHisto; 
    TH2F * effDenHisto; 
    
    evtHisto = (TH1F*) hitFile->GetObjectChecked("hGood", "TH1F");
    colvscluHisto = (TH2F*) cluFile->GetObjectChecked("hncDoubleColumn", "TH2F"); 
    cluHisto1 = colvscluHisto->ProjectionX("projX"); 
    nDoubleColumnsHisto = colvscluHisto->ProjectionY("projY"); 
    effNumHisto = (TH2F*) effFile->GetObjectChecked("h2_num", "TH2F"); 
    effDenHisto = (TH2F*) effFile->GetObjectChecked("h2_den", "TH2F"); 
      
    if (!effNumHisto) {
      cout << "Can't open effNumHisto!" << endl; 
      continue; 
    }
    if (!effDenHisto) {
      cout << "Can't open effDenHisto!" << endl; 
	continue; 
    }
    
    TH1F* hEmpty = (TH1F*) hitFile->GetObjectChecked("hEmpty", "TH1F");
   
    if (!(pixHisto&&cluHisto1&&evtHisto)) {
	std::cout << Form("One (or more) of the required objects for spill %d is missing from its file", spill) << std::endl;
	continue;
    }
    int nSlices, iBin, iDC;

    nSlices = minValue(pixHisto->GetNbinsX(), cluHisto1->GetNbinsX());
    
    nSlices = minValue(nSlices, evtHisto->GetNbinsX());
    int nDCs ; // Number of Double Columns 
    nDCs = nDoubleColumnsHisto->GetNbinsX(); 

    for (slice=0; slice<nSlices; ++slice) {
      iBin = slice+1;
      nHitPixels = pixHisto->GetBinContent(iBin);

      nEvents = evtHisto->GetBinContent(iBin);
      nEmptyEvents = hEmpty->GetBinContent(iBin);
      for (DC=0; DC < nDCs; ++DC) {
	iDC = DC + 1; 
	nClusters = colvscluHisto->GetBinContent(iBin, iDC);
	eff_num = effNumHisto->GetBinContent(iBin, iDC); 
	eff_den = effDenHisto->GetBinContent(iBin, iDC); 
	summary->Fill();

      if ( label == "v4") {
	for (DC=0; DC < nDCs; ++DC) {
	  iDC = DC + 1; 
	  // nClusters = colvscluHisto->GetBinWithContent2(0, iBin, iDC);
	  // Int_t bin = colvscluHisto->GetBin(iBin, iDC);
	  // eff_num = effNumHisto->GetBinContent(bin); 
	  // eff_den = effDenHisto->GetBinContent(bin); 

	  nClusters = colvscluHisto->GetBinContent(iBin, iDC);
	  eff_num = effNumHisto->GetBinContent(iBin, iDC); 
	  eff_den = effDenHisto->GetBinContent(iBin, iDC); 	  

	  // cout << "slice = " << slice << " , DC = " << DC 
	  //      << " , nClusters = " << nClusters
	  //      << " , eff_num = " << eff_num 
	  //      << " , eff_den = " << eff_den 
	  //      << endl ; 
	  
	  summary->Fill();
	}
      }
      }
    }
    
    pixFile->Close();
    cluFile->Close();
    evtFile->Close();
    effFile->Close();
    hitFile->Close();
  
  } // Loop over logbook spills

  summary->SaveAs(outfile);
}




#ifndef __CINT__ 
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
  cerr << "Usage: ./summary label \n" << endl; 
}

int main(int argc, char** argv) {
  if (argc < 2) {
    print_usage() ;  
    return -1; 
  }

  TString label = argv[1]; 
  if (label != "roc2") {
    cerr << "Please use ./summary roc2 !" << endl; 
    return -1; 
  }

  summary(label); 
  
  gSystem->Exit(0);

  return 0 ;
}

#endif

