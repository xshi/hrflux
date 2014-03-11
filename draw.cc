// Author: Xin Shi <Xin.Shi@cern.ch>
// Based on Stefano Mersi's helper.cpp 
// Created: [2013-03-28 Thu 08:39] 


#include <TObject.h>
#include <TColor.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <string>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <TF1.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <math.h>
//#include <Palette.h>

#include <stdlib.h>
using namespace std; 

const Long64_t unknownLines = -1; // Must be negative
const double sqrt_2_pi = 2.50662827463100024;
const double limitProbability = 1e-12;

/*****************************************/
/*                                       */
/*          Poisson magic trick          */
/*                                       */
/*****************************************/
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

  // gStyle->SetPadTopMargin(PadTopMargin);  
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

TPaveText * get_FitParameters(Double_t x1, Double_t y1, Double_t x2, Double_t y2, 
			      TF1 *f1, int npar){
  TPaveText *pt = new TPaveText(x1, y1, x2, y2);
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

string get_cut(TString label){
  string cut = ""; 
  // if (ver == 2) {
  if (label == "v2" 
      or label == "v3" )
    {
      cut = "slice>2&&slice<18&&nClusters>nEvents&&flux_roc<800";
    }
  
  if (label == "v4" ) {
    // cut = "slice>2&&slice<18&&nClusters>nEvents";
    cut = "slice>2&&slice<18";
  }
  // cout << "cut = " << cut << endl; 
  return cut; 
}

// This gives the real lambda of a poissionian whose zero-th bin was
// omitted from the counts obtaining an average of m
double invertProb(double m) {
  double lambda;
  if (m<=1) return 0;

  if (m<1.5) {
    lambda = 2* (1-1/m);
    if (m<1+1e-6) return lambda;
  } else {
    lambda = m;
  }

  double x;
  for (int i=0; i<5; ++i) {
    x = (1-exp(-lambda) - lambda/m)/(1/m-exp(-lambda));
    lambda = lambda+x;
  }

  return lambda;
}

// This gives the measured average of a poissionian with real average
// x, but whose zero-th bin was omitted from the counts
double directProb(double x) {
  return  x/(1-exp(-x));
}

Color_t PaletteColor(const unsigned int& plotIndex) {
  std::string colorCode;
  
  if (plotIndex==0) colorCode = "#000000";
  else {
    int nColor=(plotIndex-1) % 12;
    switch (nColor) {
    case 0 :
      colorCode="#004586";
      break;
    case 1 :
      colorCode="#FF420E";
      break;
    case 2 :
      colorCode="#FFD320";
      break;
    case 3 :
      colorCode="#579D1C";
      break;
    case 4 :
      colorCode="#7E0021";
      break;
    case 5 :
      colorCode="#83CAFF";
      break;
    case 6 :
      colorCode="#314004";
      break;
    case 7 :
      colorCode="#AECF00";
      break;
    case 8 :
      colorCode="#4B1F6F";
      break;
    case 9 :
      colorCode="#FF950E";
      break;
    case 10 :
      colorCode="#C5000B";
      break;
    case 11 :
      colorCode="#0084D1";
      break;
    default :
      std::cerr << "ERROR: in Vizard::getNiceColor() n%12 is not an int between 0 and 11! This should not happen." << std::endl;
      colorCode="#000000";
      break;
    }
  }
  return TColor::GetColor(colorCode.c_str());
}



bool createVector(TTree* summary, const std::string& xval,
                  const std::string& selection, Long64_t& expectedLines, Double_t* &result,
                  const Double_t* multiplier = NULL) {
  summary->Draw(xval.c_str(), selection.c_str(), "goff");
  Double_t* values = summary->GetV1();
  Long64_t nLines = summary->GetSelectedRows();
  if (expectedLines==unknownLines) expectedLines=nLines;
  if (nLines!=expectedLines) {
    return false;
  }

  result = (Double_t*)malloc(sizeof(Double_t)*nLines);
  if (xval=="") for (int i=0; i<nLines; ++i) result[i] = 0;
  else if (multiplier==NULL) for (int i=0; i<nLines; ++i) result[i] = values[i];
  else for (int i=0; i<nLines; ++i) result[i] = values[i] * multiplier[i];
  return true;
}

TGraphAsymmErrors* createGraph(TTree* summary,
                               const std::string& xval, const std::string& yval,
                               const std::string& selection, 
                               const std::string& xerr_lo, const std::string& xerr_hi,
                               const std::string& yerr_lo, const std::string& yerr_hi,
                               bool xerrRelative = false, bool yerrRelative = false) {

  summary->SetEstimate(summary->GetEntries());
  
  Double_t *xx, *yy, *exl, *exh, *eyl, *eyh;
  bool okxx, okyy, okexl, okexh, okeyl, okeyh;
  Long64_t nLines = unknownLines;
  // X values
  okxx = createVector(summary, xval, selection, nLines, xx); if (!okxx) return NULL;
  // Y values
  okyy = createVector(summary, yval, selection, nLines, yy); if (!okyy) return NULL;
  // X error (low)
  if (xerrRelative==true) okexl = createVector(summary, xerr_lo, selection, nLines, exl, xx);
  else okexl = createVector(summary, xerr_lo, selection, nLines, exl);
  if (!okexl) return NULL;
  // X error (high, if needed)
  if (xerr_hi=="") {
    exh = exl;
    okexh = true;
  } else {
    if (xerrRelative==true) okexh = createVector(summary, xerr_hi, selection, nLines, exh, xx);
    else okexh = createVector(summary, xerr_hi, selection, nLines, exh);
  }
  if (!okexh) return NULL;
  // Y error (low)
  if (yerrRelative==true) okeyl = createVector(summary, yerr_lo, selection, nLines, eyl, yy);
  else okeyl = createVector(summary, yerr_lo, selection, nLines, eyl);
  if (!okeyl) return NULL;
  // Y error (high, if needed)
  if (yerr_hi=="") {
    eyh = eyl;
    okeyh = true;
  } else {
    if (yerrRelative==true) okeyh = createVector(summary, yerr_hi, selection, nLines, eyh, yy);
    else okeyh = createVector(summary, yerr_hi, selection, nLines, eyh);
  }
  if (!okeyh) return NULL;

  TGraphAsymmErrors* result = new TGraphAsymmErrors(nLines, xx, yy, exl, exh, eyl, eyh);
  return result;
}

double gaussian(double x_over_sigma) {
  double nn = x_over_sigma*x_over_sigma;
  nn = exp(-0.5*nn);
  return nn;//*sqrt_2_pi;
}

TH2D* probabilityDensity(TGraphAsymmErrors* myGraph,
                         const std::string name, const std::string title,
                         Int_t xBins, Double_t xLow, Double_t xHigh,
                         Int_t yBins, Double_t yLow, Double_t yHigh,
                         Int_t nPointsScan=100, Double_t nSigmas = 3, 
                         bool flatx = false, bool flaty = false) {

  TH2D* result = new TH2D(name.c_str(), title.c_str(),
                          xBins, xLow, xHigh,
                          yBins, yLow, yHigh);
  
  
  Double_t *xx, *yy, *exl, *exh, *eyl, *eyh;
  xx = myGraph->GetX();
  yy = myGraph->GetY();
  exl = myGraph->GetEXlow();
  exh = myGraph->GetEXhigh();
  eyl = myGraph->GetEYlow();
  eyh = myGraph->GetEYhigh();
  Double_t px, py, weight;
  weight = 1/double(2*nPointsScan+1);
  weight *= weight ;
  double wx, wy, ex, ey;
  double nps = double(nPointsScan);
  TProfile2D* pointResult[myGraph->GetN()];
  for (Int_t iPoint=0; iPoint<myGraph->GetN(); ++iPoint) {
    pointResult[iPoint] = new TProfile2D(Form("pointImprob%d", iPoint), "Improb",
                                         xBins, xLow, xHigh,
                                         yBins, yLow, yHigh);
    for (int j=-nPointsScan; j<=nPointsScan; ++j) {
      if (j<=0) ex = exl[iPoint]; else ex = exh[iPoint];
      px = j * ex * nSigmas/nps;
      if (flatx) wx = 1; else wx = gaussian(j/nps*nSigmas);
      px += xx[iPoint];
      for (int k=-nPointsScan; k<=nPointsScan; ++k) {
        if (k<=0) ey = eyl[iPoint]; else ey = eyh[iPoint];
        py = k * ey * nSigmas/nps;
        py += yy[iPoint];
        if (flaty) wy = 1; else wy = gaussian(k/nps*nSigmas);
        // We can fill px, py
        pointResult[iPoint]->Fill(px, py, wx*wy*weight);
      }
    }
  }


  double value;
  for (int i=1; i<=result->GetNbinsX(); ++i) {
    px = result->GetXaxis()->GetBinCenter(i);
    for (int j=1; j<=result->GetNbinsY(); ++j) {
      py = result->GetYaxis()->GetBinCenter(j);
      weight = 1;
      for (int k=0; k<myGraph->GetN(); ++k) {
        value = pointResult[k]->GetBinContent(i,j);
        if (value>limitProbability) {
          weight *= 1-value;
        }
      }
      value = 1-weight;
      if (value>limitProbability) result->SetBinContent(i, j, value);
    } 
  }

  for (int k=0; k<myGraph->GetN(); ++k) delete pointResult[k];

  return result;
}

std::string golden0 ="WBC==100&&HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800";
std::string golden1 ="WBC==100&&HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0";
std::string golden2 = "WBC==100&&HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>1";
std::string golden3 = "WBC==100&&HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>1&&flux_cm>0&&slice<20";
std::string golden3WBCFree = "HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>1&&flux_cm>0&&slice<20";
std::string golden4 = "WBC==100&&HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>2&&flux_cm>0&&slice<18";
std::string golden4WBCFree = "HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>2&&flux_cm>0&&slice<18";
std::string psFree = "WBC==100&&HV==-120&&vthrcomp==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>1&&flux_cm>0";
std::string vthrcompFree = "WBC==100&&HV==-120&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>1&&flux_cm>0";
std::string ps_d = "WBC==100&&HV==-120&&vthrcomp==0&&PS==%d&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>1&&flux_cm>0";
std::string wbcFree = "HV==-120&&vthrcomp==0&&PS==0&&x_mm==-2&&y_pos==2800&&nClusters!=0&&nHitPixels!=0&&nEvents!=0&&slice>1&&flux_cm>0";
std::string minimumSelection =  "slice>2&&slice<18&&nClusters>nEvents";

double rocArea = 0.648; 

TGraphAsymmErrors* gg;
TH2D* pippo;
TProfile* ll;
TCanvas* myCanvas;
TTree* summary;

TGraph* comparisonEfficiency = NULL;

void makeEffVsWBC() {
  gg = createGraph(summary, "WBC", "efficiency*100", wbcFree, "", "", "", "");
  gg->SetMarkerStyle(7);
  gg->SetMarkerColor(PaletteColor(1));
  gg->SetTitle("Efficiency;WBC [CK];Efficiency [%]");
  myCanvas = new TCanvas();
  gg->Draw("ap");
}

void fluxTimes() {
  gg = createGraph(summary, "time", "flux_roc", golden3, "0.5", "", "flux_relerr", "", false, true);
  TH2D* frame = new TH2D("frame", "Flux on ROC;Time [s];Flux [MHz]", 10, 0, 12, 10, 0, 900);
  frame->Draw();
  gg->Draw("p same");
  // pippo = probabilityDensity(gg, "pippo",  "Flux on ROC;Time [s];Flux [MHz]", 100, 0, 12, 100, 0, 900, 40, 3);
  // pippo->Draw("colz");
}

void efficiencyVsHitsPerSlit() {
  std::string slitSelection[2];
  TF1* myPol[2];
  TF1* myPolFlux[2];
  TGraphAsymmErrors* plots[2];
  TGraphAsymmErrors* plotsFlux[2];
  int slits[] = {20, 40};
  int plotColors[] = {2,0};
  //int markerStyles[] = {8, 7};

  double scaleFactor= 1.22510;

  myCanvas = new TCanvas();
  myCanvas->Divide(2,1);

  TH2D* frame1 = new TH2D("frame1", "Efficiency;Cluster mutiplicity;Efficiency [%]", 10, 0, 10, 10, 0.5, 100);
  TH2D* frame2 = new TH2D("frame2", "Efficiency;Hit-derived flux [MHz/cm^{2}];Efficiency [%]", 10, 0, 350/rocArea*scaleFactor, 10, 0.5, 100);
  myCanvas->cd(1);
  frame1->Draw();
  myCanvas->cd(2);
  frame2->Draw();  
  
  std::string mySelection = minimumSelection + "&&magslits_p==%d" + "&&efficiency>0";
  double maxHits = 5; // Maximum x value for the linear fit in nHits

  for (int i=0; i<2; ++i) {
    myPol[i]= new TF1("myPol1", "pol1", 0, 5), myPol[i]->SetLineColor(kBlue);
    myPolFlux[i]= new TF1("myPolFlux", "pol1", 0, 250), myPolFlux[i]->SetLineColor(kBlue);
    slitSelection[i] = Form(mySelection.c_str(), slits[i]);

    myCanvas->cd(1);
    plots[i] = createGraph(summary, "invertProb(nClusters/nEvents)", "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    plots[i]->SetLineColor(PaletteColor(plotColors[i]));
    plots[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plots[i]->Draw("p same");
    plots[i]->Fit(myPol[i], "S", "", 0, maxHits);

    myCanvas->cd(2);

    plotsFlux[i] = createGraph(summary, Form("invertProb(nClusters/nEvents)*40/%f*%f", rocArea, scaleFactor), "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    plotsFlux[i]->SetLineColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->Draw("p same");
    plotsFlux[i]->Fit(myPolFlux[i], "S", "", 0, maxHits*40/rocArea*scaleFactor);
    // summary->SetMarkerColor(PaletteColor(plotColors[i]));
    // summary->SetMarkerStyle(markerStyles[i]);
    // summary->Draw("efficiency:nClusters/nEvents*40", slitSelection[i].c_str(), "same p");
  }
  for (int i=0; i<2; ++i) {
   // Print results here
  }
comparisonEfficiency->Draw("same l");
}

/*
  // How to use the above function:
  .L helper.cpp+
  initialize();
  efficiencyVsHitsPerSlit();
  myCanvas->cd(2);
  comparisonEfficiency->SetLineColor(kGreen);
  comparisonEfficiency->Draw("same l");
*/

void makeFluxHistogram() {
  TH2D* frame = new TH2D("frame", "Flux on ROC;Flux [MHz];Counts", 100, 0, 600, 10, 0, 140);
  frame->Draw();
  summary->SetLineWidth(2);
  summary->SetLineColor(PaletteColor(2));
  summary->Draw("flux_roc", golden3.c_str(), "same");
  summary->SetLineColor(PaletteColor(1));
  summary->Draw("(nClusters/nEvents)*40", golden3.c_str(), "same");
}

void makeEffvsFlux() {
  gg = createGraph(summary, "flux_roc",  "efficiency*100", golden3, "flux_relerr", "", "eff_lo*100", "eff_hi*100", true, false);
  gg->SetTitle("Efficiency;Flux on ROC [MHz];Efficiency [%]");
  gg->Draw("ap");
}


void makeEffvsFluxSliced() {
  std::string golden3Sliced[20];
  int i;
  for (i=0; i<20; ++i) { golden3Sliced[i]=Form("%s&&slice==%d", golden3.c_str(), i+1);}
  myCanvas = new TCanvas();
  myCanvas->Divide(5,4);
  TH2D* frame;

  TCanvas* otherCanvas = new TCanvas();
  frame = new TH2D("frameall","Efficiency in all slices;Flux on ROC [MHz];Efficiency [%%]" , 100, 0, 900, 10, 0, 100);
  frame->Draw();

  for (i=0; i<20; ++i) {
    myCanvas->cd(i+1);
    gg = createGraph(summary, "flux_roc", "efficiency*100", golden3Sliced[i], "flux_relerr", "", "eff_lo*100", "eff_hi*100", true, false);
    gg->SetTitle(Form("Efficiency in slice %d;Flux on ROC [MHz];Efficiency [%%]", i+1));
    gg->SetLineColor(PaletteColor(i-1));
    if (gg->GetN()>0) {
      //frame = new TH2D(Form("frame%d", i), gg->GetTitle(), 100, 0, 900, 10, 0, 100);
      //frame->Draw();
      gg->Draw("ap");
      
      otherCanvas->cd();
      gg->Draw("p same");
    }
  }  
}

void makeHitsVsFlux() {
  gg = createGraph(summary, "flux_roc", "nHitPixels/nEvents", golden4, "flux_relerr", "", "1/sqrt(nHitPixels)", "", true, true);
  gg->SetTitle("N Hits;Flux on ROC [MHz];Average N hits");
  gg->Draw("ap");
  TF1* myPol1 = new TF1("myPol1", "pol1", 100, 400);
  gg->Fit(myPol1, "S", "", 0, 300);
  myPol1->SetLineColor(kRed);
  myPol1->Draw("same");
}

void makeClustersVsFluxWBC() {
  std::vector<int> WBCList;
  std::vector<int>::iterator it;
  std::string selection;

  WBCList.push_back(100);
  WBCList.push_back(20);
  WBCList.push_back(30);
  WBCList.push_back(40);
  WBCList.push_back(50);
  WBCList.push_back(60);
  WBCList.push_back(70);
  WBCList.push_back(80);
  WBCList.push_back(110);
  WBCList.push_back(120);

  myCanvas = new TCanvas();
  myCanvas->Divide(4,5);

  int iColor;
  std::string plotOpt;
  TH2D* frame;
  for (int i=1; i<=20; ++i) {
    myCanvas->cd(i);
    plotOpt = "ap";
    iColor=0;
    frame = new TH2D("frameall", "N Clusters;Flux on ROC [MHz];Average N hits", 100, 0, 400, 10, 0, 10);
    frame->Draw();
    plotOpt = "p same";
    for (it=WBCList.begin(); it!=WBCList.end(); ++it) {
      selection = golden3WBCFree + Form("&&WBC==%d&&slice==%d", *it, i);
      gg = createGraph(summary, "flux_roc", "nClusters/nEvents", selection, "flux_relerr", "", "1/sqrt(nClusters)", "", true, true);
      //gg->SetTitle("N Clusters;Flux on ROC [MHz];Average N hits");
      gg->SetLineColor(PaletteColor(iColor));
      if (gg->GetN()) {
        gg->Draw(plotOpt.c_str());
        plotOpt = "p same";
      }
      iColor++;
    }
  }
}

void makeClustersVsFlux() {
  gg = createGraph(summary, "flux_roc", "nClusters/nEvents", golden4, "flux_relerr", "", "1/sqrt(nHitPixels)", "", true, true);
  gg->SetTitle("N Clusters;Flux on ROC [MHz];Average N hits");
  gg->Draw("ap");  
  TF1* myPol1 = new TF1("myPol1", "pol1", 100, 400);
  gg->Fit(myPol1, "S", "", 0, 300);
  myPol1->SetLineColor(kRed);
  myPol1->Draw("same");
}

/*
   // Uncool
void makeClustersVsFluxSlits() {
  TH2D* frame = new TH2D("frame", "Cluster multiplicity vs flux;Flux on ROC [MHz];nClusters/nEvents", 10, 0, 400, 10, 0, 10);
  frame->Draw();
  summary->SetMarkerColor(kRed);
  summary->Draw("nClusters/nEvents:flux_roc", "flux_roc>0&&flux_roc<400&&magslits_p==20&&slice>2&&nClusters>nEvents&&slice<18", "p same");
  summary->SetMarkerColor(kBlack);
  summary->Draw("nClusters/nEvents:flux_roc", "flux_roc>0&&flux_roc<400&&magslits_p==40&&slice>2&&nClusters>nEvents&&slice<18", "p same");
}
*/
// COOL
void makeClustersVsFluxSlits() {
  std::string slitSelection[2];
  TF1* myPol[2];
  TGraphAsymmErrors* plots[2];
  int slits[] = {20, 40};
  int plotColors[] = {2,0};

  TH2D* frame = new TH2D("frame", "Cluster multiplicity vs flux;Flux on ROC [MHz];nClusters/nEvents", 10, 0, 400, 10, 0, 10);
  frame->Draw();
  for (int i=0; i<2; ++i) {
    myPol[i]= new TF1("myPol1", "pol1", 100, 200), myPol[i]->SetLineColor(kBlue);
    slitSelection[i] = Form("flux_roc>0&&flux_roc<400&&magslits_p==%d&&slice>2&&nClusters>nEvents&&slice<18", slits[i]);
    plots[i] = createGraph(summary, "flux_roc", "nClusters/nEvents", slitSelection[i], "flux_relerr", "", "1/sqrt(nHitPixels)", "", true, true);
    plots[i]->SetLineColor(PaletteColor(plotColors[i]));
    plots[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plots[i]->Draw("p same");
    plots[i]->Fit(myPol[0], "S", "", 0, 200);
  }
  for (int i=0; i<2; ++i) {
   // Print results here
  }
}


void makeClustersVsFluxCorrected() {
  gg = createGraph(summary, "flux_roc", "(nClusters/nEvents)", golden4, "flux_relerr", "", "1/sqrt(nHitPixels)", "", true, true);
  gg->SetTitle("N Clusters;Flux on ROC [MHz];Average N hits");
  gg->Draw("ap");  
  TF1* myPol1 = new TF1("myPol1", "pol1", 100, 300);
  gg->Fit(myPol1, "S", "", 0, 300);
  myPol1->SetLineColor(kRed);
  myPol1->Draw("same");

  double offset = myPol1->GetParameter(0);
  //double slope = myPol1->GetParameter(1);

  myCanvas = new TCanvas();
  TH2D* frame = new TH2D("frameall", "AAAAAAAAA", 10, -1, 1, 10, 0, 100);
  frame->Draw();
  
  // std::vector<int> WBCList;
  // std::vector<int>::iterator it;
  // WBCList.push_back(100);
  // WBCList.push_back(20);
  // WBCList.push_back(30);
  // WBCList.push_back(40);
  // WBCList.push_back(50);
  // WBCList.push_back(60);
  // WBCList.push_back(70);
  // WBCList.push_back(80);
  // WBCList.push_back(110);
  // WBCList.push_back(120);

  std::string plotOpt = "colz";
  std::string selection;
  //  int iColor=0;
  //  for (it=WBCList.begin(); it!=WBCList.end(); ++it) {
  // selection = golden4WBCFree + Form("&&WBC==%d", *it) + Form("&&(nClusters/nEvents)>%f", offset);

  selection = golden4WBCFree + Form("&&(nClusters/nEvents)>%f", offset);
  //summary->SetLineColor(PaletteColor(iColor));
  summary->Draw(Form("((nClusters/nEvents)-%f)/flux_roc:WBC", offset), selection.c_str(), plotOpt.c_str());
  plotOpt = "same";
  //iColor++;
  
    //  }

}

void makeCorrectedEfficiency() {
  gg = createGraph(summary, "invertProb(nClusters/nEvents)*40", "efficiency*100", golden3, "1/sqrt(nClusters)", "", "eff_lo*100", "eff_hi*100", true, false);
  gg->SetTitle("Efficiency;Hits-derived flux on ROC [MHz];Efficiency [%]");
  gg->Draw("ap");
}

void initialize(TString label) {
  TString infile = Form("/raid1/w/xshi/summary/summary_%s.root", label.Data()); 
  TFile *_file0 = TFile::Open(infile);
  summary = (TTree*) _file0->GetObjectChecked("summary", "TTree");

  if (!comparisonEfficiency) {
    comparisonEfficiency = new TGraph();
    comparisonEfficiency->SetName("comparisonEfficiency");
    int iPoint=0;
    comparisonEfficiency->SetPoint(iPoint++, 15, 99.5);
    comparisonEfficiency->SetPoint(iPoint++, 30, 99.5);
    comparisonEfficiency->SetPoint(iPoint++, 50, 99.4);
    comparisonEfficiency->SetPoint(iPoint++, 63, 99.3);
    comparisonEfficiency->SetPoint(iPoint++, 116, 98.8);
    comparisonEfficiency->SetPoint(iPoint++, 142, 98.4);
    comparisonEfficiency->SetPoint(iPoint++, 267, 94.4);
    comparisonEfficiency->SetLineWidth(2);
    comparisonEfficiency->SetLineColor(kBlue);
    comparisonEfficiency->Draw("same");
  }

}

void drawFluxHistogram(TString label, TString figfile) {
  // summary->Draw("flux_roc", "flux_roc<500"); 
  // TCanvas *c = new TCanvas("c", "c", 600, 600);

  TH2D* frame = new TH2D("frame", "Flux on ROC;Flux [MHz];Counts", 100, 0, 600, 10, 0, 140);
  frame->Draw();
  summary->SetLineWidth(2);
  summary->SetLineColor(PaletteColor(2));
  // summary->Draw("flux_roc", "flux_roc<500", "same"); 
  // summary->Draw("flux_roc", golden0, "same"); 


  // string cut = "flux_roc<500"; 
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 

  summary->Draw("flux_roc", cut.c_str(), "same"); 


  // TString figpath = "/afs/cern.ch/user/x/xshi/www/pxl/fig"; 
  // TString figfile;
  // figfile.Form("%s/EffvsFluxSliced_v%d.pdf", figpath.Data(), ver); 
  // figfile.Form("%s/FluxHistogram_v%d.pdf", figpath.Data(), ver); 
  myCanvas->Print(figfile);
  // delete c; 
  // summary->Draw("flux_roc", golden3.c_str(), "same");
  // summary->SetLineColor(PaletteColor(1));
  // summary->Draw("(nClusters/nEvents)*40", golden3.c_str(), "same");
}

void drawEffvsFlux(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 
  gg = createGraph(summary, "flux_roc",  "efficiency*100", cut, "flux_relerr", "", "eff_lo*100", "eff_hi*100", true, false);
  gg->SetTitle("Efficiency;Flux on ROC [MHz];Efficiency [%]");


  // TCanvas *c = new TCanvas("c", "c", 600, 600);
  
  gg->Draw("ap");
  
  // TString figpath = "/afs/cern.ch/user/x/xshi/www/pxl/fig"; 
  // TString figfile;
  // figfile.Form("%s/EffvsFlux_v%d.pdf", figpath.Data(), ver); 
  myCanvas->Print(figfile);
  // delete c; 
 
}

void drawEffvsFluxSliced(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 
  std::string cutSliced[20];
  int i;
  for (i=0; i<20; ++i) { cutSliced[i]=Form("%s&&slice==%d", cut.c_str(), i+1);}
  TH2D* frame;
  frame = new TH2D("frameall","Efficiency in all slices;Flux on ROC [MHz];Efficiency [%]" , 100, 0, 600, 10, 0, 100);
  frame->Draw();

  for (i=0; i<20; ++i) {
    gg = createGraph(summary, "flux_roc", "efficiency*100", cutSliced[i], "flux_relerr", "", "eff_lo*100", "eff_hi*100", true, false);
    gg->SetTitle(Form("Efficiency in slice %d;Flux on ROC [MHz];Efficiency [%%]", i+1));
    gg->SetLineColor(PaletteColor(i-1));
    if (gg->GetN()>0) {
      gg->Draw("p same");
   }
  }
  myCanvas->Print(figfile);
}


void drawClustersVsFluxSlits(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 
  std::string slitSelection[2];
  TF1* myPol[2];
  TGraphAsymmErrors* plots[2];
  int slits[] = {20, 40};
  int plotColors[] = {2,0};

  TH2D* frame = new TH2D("frame", "Cluster multiplicity vs flux;Flux on ROC [MHz];nClusters/nEvents", 10, 0, 400, 10, 0, 10);
  frame->Draw();
  for (int i=0; i<2; ++i) {
    myPol[i]= new TF1("myPol1", "pol1", 100, 200), myPol[i]->SetLineColor(kBlue);
    // slitSelection[i] = Form("flux_roc>0&&flux_roc<400&&magslits_p==%d&&slice>2&&nClusters>nEvents&&slice<18", slits[i]);
    slitSelection[i] = Form("%s&&slits==%d", cut.c_str(), slits[i]);

   if (label == "v2") 
    plots[i] = createGraph(summary, "flux_roc", "nClusters/nEvents", slitSelection[i], "flux_relerr", "", "1/sqrt(nHitPixels)", "", true, true);

   else if (label == "v3") 
     plots[i] = createGraph(summary, "flux_roc", "nClusters/(nEvents+nEmptyEvents)", slitSelection[i], "flux_relerr", "", "1/sqrt(nHitPixels)", "", true, true);

   else return; 
   
 
    plots[i]->SetLineColor(PaletteColor(plotColors[i]));
    plots[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plots[i]->Draw("p same");
    plots[i]->Fit(myPol[0], "S", "", 0, 200);
  }

  myCanvas->Print(figfile);

  for (int i=0; i<2; ++i) {
   // Print results here
  }
}


void drawEffvsHitsPerSlit(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 

  std::string slitSelection[2];
  TF1* myPol[2];
  TF1* myPolFlux[2];
  // TGraphAsymmErrors* plots[2];
  TGraphAsymmErrors* plotsFlux[2];
  int slits[] = {20, 40};
  int plotColors[] = {2,0};
  //int markerStyles[] = {8, 7};

  double scaleFactor= 1.22510;

  // myCanvas = new TCanvas();
  // myCanvas->Divide(2,1);

  // TH2D* frame1 = new TH2D("frame1", "Efficiency;Cluster mutiplicity;Efficiency [%]", 10, 0, 10, 10, 0.5, 100);
  TH2D* frame2 = new TH2D("frame2", "Efficiency;Hit-derived flux [MHz/cm^{2}];Efficiency [%]", 10, 0, 350/rocArea*scaleFactor, 10, 0.5, 100);
  // myCanvas->cd(1);
  // frame1->Draw();
  // myCanvas->cd(2);
  frame2->Draw();  
  
  // std::string mySelection = minimumSelection + "&&magslits_p==%d" + "&&efficiency>0";
  // std::string mySelection = minimumSelection; 
  double maxHits = 5; // Maximum x value for the linear fit in nHits

  for (int i=0; i<2; ++i) {
    myPol[i]= new TF1("myPol1", "pol1", 0, 5), myPol[i]->SetLineColor(kBlue);
    myPolFlux[i]= new TF1("myPolFlux", "pol1", 0, 250), myPolFlux[i]->SetLineColor(kBlue);
    // slitSelection[i] = Form(mySelection.c_str(), slits[i]);
    slitSelection[i] = Form(cut.c_str(), slits[i]);

    // myCanvas->cd(1);
    // plots[i] = createGraph(summary, "invertProb(nClusters/nEvents)", "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    // plots[i]->SetLineColor(PaletteColor(plotColors[i]));
    // plots[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    // plots[i]->Draw("p same");
    // plots[i]->Fit(myPol[i], "S", "", 0, maxHits);

    // myCanvas->cd(2);

    plotsFlux[i] = createGraph(summary, Form("invertProb(nClusters/nEvents)*40/%f*%f", rocArea, scaleFactor), "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    plotsFlux[i]->SetLineColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->Draw("p same");
    plotsFlux[i]->Fit(myPolFlux[i], "S", "", 0, maxHits*40/rocArea*scaleFactor);
    // summary->SetMarkerColor(PaletteColor(plotColors[i]));
    // summary->SetMarkerStyle(markerStyles[i]);
    // summary->Draw("efficiency:nClusters/nEvents*40", slitSelection[i].c_str(), "same p");
  }

  myCanvas->Print(figfile);

  for (int i=0; i<2; ++i) {
   // Print results here
  }
  // comparisonEfficiency->Draw("same l");
}



void drawHitsVsFlux(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 

  gg = createGraph(summary, "flux_roc", "nHitPixels/nEvents", cut.c_str(), "flux_relerr", "", "1/sqrt(nHitPixels)", "", true, true);
  gg->SetTitle("N Hits;Flux on ROC [MHz];Average N hits");
  gg->Draw("ap");
  TF1* myPol1 = new TF1("myPol1", "pol1", 100, 400);
  gg->Fit(myPol1, "S", "", 0, 300);
  myPol1->SetLineColor(kRed);
  myPol1->Draw("same");


  TPaveText *pt = new TPaveText(300, 2, 400, 4);

  pt->SetBorderSize(0);

  pt->SetFillColor(0); 
  Char_t message[80];
  sprintf(message, "slope = %.2f ns", myPol1->GetParameter(1)*1000);
  pt->AddText(message);  
  pt->Draw();
  myCanvas->Print(figfile);

}

void drawEffvsClusterMultiplicity(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 

  std::string slitSelection[2];
  TF1* myPol[2];
  TF1* myPolFlux[2];
  TGraphAsymmErrors* plots[2];
  // TGraphAsymmErrors* plotsFlux[2];
  int slits[] = {20, 40};
  int plotColors[] = {2,0};
  //int markerStyles[] = {8, 7};

  // double scaleFactor= 1.22510;

  // myCanvas->cd();
  // myCanvas = new TCanvas("c", "c", 1200, 600);
  // myCanvas = new TCanvas("c", "c", 600, 600);
  // myCanvas->Divide(2,1);

  TH2D* frame1 = new TH2D("frame1", "Efficiency;Cluster mutiplicity;Efficiency [%]", 10, 0, 10, 10, 0.5, 100);
  //   TH2D* frame2 = new TH2D("frame2", "Efficiency;Hit-derived flux [MHz/cm^{2}];Efficiency [%]", 10, 0, 350/rocArea*scaleFactor, 10, 0.5, 100);
  // myCanvas->cd(1);
  frame1->Draw();
  //  myCanvas->cd(2);
  // frame2->Draw();  
  
  std::string mySelection = minimumSelection + "&&magslits_p==%d" + "&&efficiency>0";
  double maxHits = 5; // Maximum x value for the linear fit in nHits

  for (int i=0; i<2; ++i) {
    myPol[i]= new TF1("myPol1", "pol1", 0, 5), myPol[i]->SetLineColor(kBlue);
    myPolFlux[i]= new TF1("myPolFlux", "pol1", 0, 250), myPolFlux[i]->SetLineColor(kBlue);
    slitSelection[i] = Form(mySelection.c_str(), slits[i]);

    myCanvas->cd(1);
    // plots[i] = createGraph(summary, "invertProb(nClusters/nEvents)", "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);

    if (label == "v2") 
      plots[i] = createGraph(summary, "(nClusters/nEvents)", "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);

    else if (label == "v3") 
      plots[i] = createGraph(summary, "(nClusters/(nEvents+nEmptyEvents))", "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);

    else return; 
    
    plots[i]->SetLineColor(PaletteColor(plotColors[i]));
    plots[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plots[i]->Draw("p same");
    plots[i]->Fit(myPol[i], "S", "", 0, maxHits);

    //     myCanvas->cd(2);

    // plotsFlux[i] = createGraph(summary, Form("invertProb(nClusters/nEvents)*40/%f*%f", rocArea, scaleFactor), "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    // plotsFlux[i] = createGraph(summary, Form("(nClusters/nEvents)*40/%f*%f", rocArea, scaleFactor), "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    // plotsFlux[i]->SetLineColor(PaletteColor(plotColors[i]));
    // plotsFlux[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    // plotsFlux[i]->Draw("p same");
    // plotsFlux[i]->Fit(myPolFlux[i], "S", "", 0, maxHits*40/rocArea*scaleFactor);
    // summary->SetMarkerColor(PaletteColor(plotColors[i]));
    // summary->SetMarkerStyle(markerStyles[i]);
    // summary->Draw("efficiency:nClusters/nEvents*40", slitSelection[i].c_str(), "same p");
  }
  for (int i=0; i<2; ++i) {
   // Print results here
  }
comparisonEfficiency->Draw("same l");

 myCanvas->Print(figfile);

}

void drawEffvsHitDerivedFlux(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 

  std::string slitSelection[2];
  TF1* myPol[2];
  TF1* myPolFlux[2];
  // TGraphAsymmErrors* plots[2];
  TGraphAsymmErrors* plotsFlux[2];
  int slits[] = {20, 40};
  int plotColors[] = {2,0};
  //int markerStyles[] = {8, 7};

  double scaleFactor= 1.22510;

  // myCanvas->cd();
  // myCanvas = new TCanvas("c", "c", 600, 600);
  // myCanvas->Divide(2,1);

  // TH2D* frame1 = new TH2D("frame1", "Efficiency;Cluster mutiplicity;Efficiency [%]", 10, 0, 10, 10, 0.5, 100);
  TH2D* frame2 = new TH2D("frame2", "Efficiency;Hit-derived flux [MHz/cm^{2}];Efficiency [%]", 10, 0, 350/rocArea*scaleFactor, 10, 0.5, 100);
  // myCanvas->cd(1);
  // frame1->Draw();
  // myCanvas->cd(2);
  frame2->Draw();  
  
  std::string mySelection = minimumSelection + "&&magslits_p==%d" + "&&efficiency>0";
  double maxHits = 5; // Maximum x value for the linear fit in nHits

  for (int i=0; i<2; ++i) {
    myPol[i]= new TF1("myPol1", "pol1", 0, 5), myPol[i]->SetLineColor(kBlue);
    myPolFlux[i]= new TF1("myPolFlux", "pol1", 0, 250), myPolFlux[i]->SetLineColor(kBlue);
    slitSelection[i] = Form(mySelection.c_str(), slits[i]);

    // myCanvas->cd(1);
    // plots[i] = createGraph(summary, "invertProb(nClusters/nEvents)", "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);

    // plots[i] = createGraph(summary, "(nClusters/nEvents)", "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    // plots[i]->SetLineColor(PaletteColor(plotColors[i]));
    // plots[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    // plots[i]->Draw("p same");
    // plots[i]->Fit(myPol[i], "S", "", 0, maxHits);

    // myCanvas->cd(2);

    // plotsFlux[i] = createGraph(summary, Form("invertProb(nClusters/nEvents)*40/%f*%f", rocArea, scaleFactor), "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);

   if (label == "v2") 
     plotsFlux[i] = createGraph(summary, Form("(nClusters/nEvents)*40/%f*%f", rocArea, scaleFactor), "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    
   else if (label == "v3") 
     plotsFlux[i] = createGraph(summary, Form("(nClusters/(nEvents+nEmptyEvents))*40/%f*%f", rocArea, scaleFactor), "efficiency*100",  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    
    else return; 

    plotsFlux[i]->SetLineColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->Draw("p same");
    plotsFlux[i]->Fit(myPolFlux[i], "S", "", 0, maxHits*40/rocArea*scaleFactor);
    // summary->SetMarkerColor(PaletteColor(plotColors[i]));
    // summary->SetMarkerStyle(markerStyles[i]);
    // summary->Draw("efficiency:nClusters/nEvents*40", slitSelection[i].c_str(), "same p");
  }
  for (int i=0; i<2; ++i) {
   // Print results here
  }
comparisonEfficiency->Draw("same l");

 myCanvas->Print(figfile);

}

void drawClustersVsDCs(TString label, TString figfile) {
  // string cut = get_cut(ver); 
  string cut = get_cut(label); 

  // summary->Draw("nClusters:DC", cut.c_str()); 
  summary->Draw("nClusters:DC"); 

  myCanvas->Print(figfile);

}

void drawEffvsHitDerivedFluxDC(TString label, TString figfile) {
  string cut = get_cut(label); 

 
  std::string slitSelection[2];
  TF1* myPol[2];
  TF1* myPolFlux[2];

  TGraphAsymmErrors* plotsFlux[2];
  int slits[] = {20, 40};
  int plotColors[] = {2,0};

  double scaleFactor= 1.22510;

  TH2D* frame2 = new TH2D("frame2", "Efficiency;Hit-derived flux [MHz/cm^{2}];Efficiency [%]", 10, 0, 350/rocArea*scaleFactor, 10, 0.5, 100);
  frame2->Draw();  
  
  std::string mySelection = minimumSelection + "&&magslits_p==%d" + "&&efficiency>0";
  double maxHits = 5; // Maximum x value for the linear fit in nHits

  for (int i=0; i<2; ++i) {
    myPol[i]= new TF1("myPol1", "pol1", 0, 5), myPol[i]->SetLineColor(kBlue);
    myPolFlux[i]= new TF1("myPolFlux", "pol1", 0, 250), myPolFlux[i]->SetLineColor(kBlue);
    slitSelection[i] = Form(mySelection.c_str(), slits[i]);
    
    plotsFlux[i] = createGraph(summary, Form("(nClusters/(nEvents+nEmptyEvents))*40/%f*%f", rocArea, scaleFactor), "(eff_num/eff_den)*100",  // string cut = get_cut(ver); 
  slitSelection[i],"1/sqrt(nClusters)+1/sqrt(nEvents)", "" , "eff_lo*100", "eff_hi*100", true, false);
    
    plotsFlux[i]->SetLineColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->SetMarkerColor(PaletteColor(plotColors[i]));
    plotsFlux[i]->Draw("p same");
    plotsFlux[i]->Fit(myPolFlux[i], "S", "", 0, maxHits*40/rocArea*scaleFactor);
    // summary->SetMarkerColor(PaletteColor(plotColors[i]));
    // summary->SetMarkerStyle(markerStyles[i]);
    // summary->Draw("efficiency:nClusters/nEvents*40", slitSelection[i].c_str(), "same p");
  }
  for (int i=0; i<2; ++i) {
   // Print results here
  }
comparisonEfficiency->Draw("same l");

 myCanvas->Print(figfile);

}


TGraphErrors * drawEffvsHitDerivedFluxSlits(TString label, TString figfile, 
					    int select_slits, bool canvas_print=true) {
  double AREA = 0.624; // cm^2     4160*150*100*10e-8
  double DC_AREA = AREA/26.;
  double TIME = 25 * pow(10, -9) ; // 25ns 
  
  // if (select_slits == 20) 
  //   TIME = 10.2 * pow(10, -9) ; // 10.2ns 
  
  const int NBINS = 50;

  double max_flux = 1000. ; // MHz 
  double BinWidth = max_flux/NBINS;
  double Nnum[NBINS];
  double Nden[NBINS];
  double Yaxis_[NBINS];
  double YaxisErr_[NBINS];
  double Xaxis_[NBINS];
  double XaxisErr_[NBINS];
  
  for(int i=0;i<NBINS;i++){
    Nnum[i] = 0;
    Nden[i] = 0;
    Xaxis_[i] = (i+0.5)*BinWidth;
    XaxisErr_[i] = (0.5)*BinWidth;
  }   
  
  
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

  Int_t nentries = (Int_t)summary->GetEntries();

  int nClustersROC = 0; 
  int nEventsROC = 0; 
  int nEmptyEventsROC = 0; 
  bool finishedOneROC = false; 
  
  int eff_num_ROC = 0; 
  int eff_den_ROC = 0; 

  summary->GetEntry(0); 
  int previous_slice = slice; 
  int nDC = 0;

  for (Int_t i=0;i<nentries;i++) {
    summary->GetEntry(i); 

    finishedOneROC = (slice!=previous_slice);

    if (finishedOneROC and nDC > 0 ) {
      int total_events; 
      total_events = nEvents + nEmptyEvents; 

      double nhits = nClustersROC*1.0/total_events; 
           
      double phi = nhits/(nDC*DC_AREA*TIME)/1000000. ; // MHz 
 
      int iBin = (int)(phi/BinWidth); 
      if ( iBin >=NBINS ) {
	continue; 
      }

      Nnum[iBin] += eff_num_ROC; 
      Nden[iBin] += eff_den_ROC;

      nClustersROC = 0;
      eff_num_ROC = 0;
      eff_den_ROC = 0;
      nEventsROC = 0; 
      nEmptyEventsROC = 0; 
      nDC = 0; 
      
    }
    
    if (slits != select_slits) continue; 

    if (slice < 4 || slice > 20 ) continue; 

    // if ( DC < 10/2 ) continue;
    // if ( DC > 40/2 ) continue;
	
    // if (eff_num < 20) continue;

    nDC++; 
    nClustersROC += nClusters; 
    eff_num_ROC += eff_num; 
    eff_den_ROC += eff_den; 
    nEventsROC += nEvents; 
    nEmptyEventsROC += nEmptyEvents; 

    previous_slice = slice;
  }

  cout << "Total entries: " << nentries << endl; 
   
  gStyle->SetOptStat(0);
  for(int i=0;i<NBINS;i++){
    if(Nden[i]!=0){
      Yaxis_[i] = 100.*Nnum[i]/Nden[i];
      YaxisErr_[i] = sqrt(Yaxis_[i]*(100.-Yaxis_[i])/Nden[i]);
    }else{
      Yaxis_[i] = 0;
      YaxisErr_[i] = 0;
    }
  }

  TGraphErrors *g = new TGraphErrors(NBINS, Xaxis_,Yaxis_,XaxisErr_,YaxisErr_);
  g->GetYaxis()->SetRangeUser(40,100);
  g->GetYaxis()->SetTitle("Efficiency [%]");
  g->GetXaxis()->SetTitle("Hit-derived flux [MHz/cm^{2}]");
  g->SetMarkerStyle(21);
  g->SetLineWidth(3);
  g->SetMarkerColor(kRed);
  // sprintf(buffer,"slits of +/- 20");
  // g->SetTitle("slits of #pm 20");
  g->SetTitle(Form("Whole ROC"));

  if (canvas_print){ 
    g->SetTitle(Form("Whole ROC: slits of #pm %d", select_slits));
    g->Draw("AEP"); 
    myCanvas->Print(figfile);
  }


  return g; 
}


TGraphErrors * drawEffvsHitDerivedFluxSlitsByDC(TString label, TString figfile, 
						int select_slits, bool canvas_print=true, 
						int secelct_DC=0, int color = 1) {
  
  double AREA = 0.624/26; // cm^2     4160*150*100*10e-8
  double TIME = 25 * pow(10, -9) ; // 25ns 
  
  // if (select_slits == 20) 
  //   TIME = 10.2 * pow(10, -9) ; // 10.2ns 
  
  const int NBINS = 50;

  double max_flux = 1000. ; // MHz 
  double BinWidth = max_flux/NBINS;
  double Nnum[NBINS];
  double Nden[NBINS];
  double Yaxis_[NBINS];
  double YaxisErr_[NBINS];
  double Xaxis_[NBINS];
  double XaxisErr_[NBINS];
  
  for(int i=0;i<NBINS;i++){
    Nnum[i] = 0;
    Nden[i] = 0;
    Xaxis_[i] = (i+0.5)*BinWidth;
    XaxisErr_[i] = (0.5)*BinWidth;
  }   
  
  
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

  Int_t nentries = (Int_t)summary->GetEntries();
  for (Int_t i=0;i<nentries;i++) {
    summary->GetEntry(i); 
     
    if (slits != select_slits) continue; 
     
    if (slice < 4 || slice > 20 ) continue; 

    // if ( DC < 1 ) continue;
    // if ( DC > 20 ) continue;

    if ( DC != secelct_DC ) continue;  
    
    int total_events = nEvents + nEmptyEvents; 
    double nhits = nClusters*1.0/total_events; 

    double phi = nhits/(AREA*TIME)/1000000. ; // MHz 
     
    int iBin = (int)(phi/BinWidth); 
    if ( iBin >=NBINS ) continue; 

    if (total_events == 0 ) 
      cout << "slits = " << slits << endl; 
     
    Nnum[iBin] += eff_num; 
    Nden[iBin] += eff_den; 
  }

  // cout << "Total entries: " << nentries << endl; 
   
  gStyle->SetOptStat(0);
  for(int i=0;i<NBINS;i++){
    if(Nden[i]!=0){
      Yaxis_[i] = 100.*Nnum[i]/Nden[i];
      YaxisErr_[i] = sqrt(Yaxis_[i]*(100.-Yaxis_[i])/Nden[i]);
    }else{
      Yaxis_[i] = 0;
      YaxisErr_[i] = 0;
    }
  }

  TGraphErrors *g = new TGraphErrors();
  int nPointsGraph;
  for (int i=0; i<NBINS; ++i) {
    nPointsGraph=g->GetN();
    if (Yaxis_[i]!=0) {
      g->SetPoint(nPointsGraph, Xaxis_[i], Yaxis_[i]);
      g->SetPointError(nPointsGraph, XaxisErr_[i], YaxisErr_[i]);
    }
  }

  // TGraphErrors *g = new TGraphErrors(NBINS, Xaxis_,Yaxis_,XaxisErr_,YaxisErr_);
  g->GetYaxis()->SetRangeUser(40,100);
  g->GetYaxis()->SetTitle("Efficiency [%]");
  g->GetXaxis()->SetTitle("Hit-derived flux [MHz/cm^{2}]");
  g->SetMarkerStyle(22);
  g->SetLineWidth(3);
  // g->SetMarkerColor(kBlue-2);
  g->SetMarkerColor(color);
  
  if (canvas_print){ 
    g->SetTitle(Form("Double Column 10 and 19: slits of #pm %d", select_slits));
    g->Draw("AEP"); 
    myCanvas->Print(figfile);
  }


  return g; 
}

TH1D * drawEffvsHitDerivedFluxSlitsByDCProfile(TString label, 
						   TString figfile, 
						   int select_slits,
						   bool canvas_print=true, 
						   int select_DC=0, 
						   int color = 1) {
  
  double AREA = 0.624/26; // cm^2     4160*150*100*10e-8
  double TIME = 25 * pow(10, -9) ; // 25ns 

  const int NBINS = 50;
  
  // TProfile *hprof  = new TProfile("hprof", "; Hit-derived flux [MHz/cm^{2}]; Efficiency [%]",
  // 				  NBINS, 0, 1000, 0, 100);
  TH1D *hprof  = new TH1D("hprof", "; Hit-derived flux [MHz/cm^{2}];",
			  NBINS, 0, 1000);

  // double max_flux = 1000. ; // MHz 
  // double BinWidth = max_flux/NBINS;
  // double Nnum[NBINS];
  // double Nden[NBINS];
  // double Yaxis_[NBINS];
  // double YaxisErr_[NBINS];
  // double Xaxis_[NBINS];
  // double XaxisErr_[NBINS];
  
  // for(int i=0;i<NBINS;i++){
  //   Nnum[i] = 0;
  //   Nden[i] = 0;
  //   Xaxis_[i] = (i+0.5)*BinWidth;
  //   XaxisErr_[i] = (0.5)*BinWidth;
  // }   
  
  
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

  Int_t nentries = (Int_t)summary->GetEntries();
  for (Int_t i=0;i<nentries;i++) {
    summary->GetEntry(i); 
     
    if (slits != select_slits) continue; 
     
    if (slice < 4 || slice > 20 ) continue; 

    // if ( DC < 1 ) continue;
    // if ( DC > 20 ) continue;

    if ( DC != select_DC ) continue;  
    
    int total_events = nEvents + nEmptyEvents; 
    double nhits = nClusters*1.0/total_events; 

    double phi = nhits/(AREA*TIME)/1000000. ; // MHz 
    
    cout << "phi = " << phi << ", eff_num / eff_den = " << eff_num << " / " << eff_den << endl; 
    


    // hprof->Fill(phi, float(eff_num)/eff_den); 
    hprof->Fill(phi); 

    // int iBin = (int)(phi/BinWidth); 
    // if ( iBin >=NBINS ) continue; 

    // if (total_events == 0 ) 
    //   cout << "slits = " << slits << endl; 
     
    // Nnum[iBin] += eff_num; 
    // Nden[iBin] += eff_den; 
  }

  // cout << "Total entries: " << nentries << endl; 
   
  gStyle->SetOptStat(0);
  // for(int i=0;i<NBINS;i++){
  //   if(Nden[i]!=0){
  //     Yaxis_[i] = 100.*Nnum[i]/Nden[i];
  //     YaxisErr_[i] = sqrt(Yaxis_[i]*(100.-Yaxis_[i])/Nden[i]);
  //   }else{
  //     Yaxis_[i] = 0;
  //     YaxisErr_[i] = 0;
  //   }
  // }

  // TGraphErrors *g = new TGraphErrors();
  // TProfile *g = new TProfile("hprof","", 100,0,400,0,100);
  // // int nPointsGraph;
  // for (int i=0; i<NBINS; ++i) {
  //   g->Fill(Xaxis_[i], Yaxis_[i]); 
  //   // nPointsGraph=g->GetN();
  //   // if (Yaxis_[i]!=0) {
  //   //   g->SetBinContent(i, Xaxis_[i], Yaxis_[i]);
  //   //   // g->SetPointError(nPointsGraph, XaxisErr_[i], YaxisErr_[i]);
  //   // }
  // }

  // TGraphErrors *g = new TGraphErrors(NBINS, Xaxis_,Yaxis_,XaxisErr_,YaxisErr_);
  // hprof->GetYaxis()->SetRangeUser(40,100);
  // hprof->GetYaxis()->SetTitle("Efficiency [%]");
  // hprof->GetXaxis()->SetTitle("Hit-derived flux [MHz/cm^{2}]");
  // hprof->SetMarkerStyle(22);
  // hprof->SetLineWidth(3);
  // // g->SetMarkerColor(kBlue-2);
  // hprof->SetMarkerColor(color);
  
  if (canvas_print){ 
    hprof->SetTitle(Form("Double Column: %d, slits: #pm %d", select_DC, select_slits));
    hprof->Draw("AEP"); 
    myCanvas->Print(figfile);
  }


  return hprof; 
}


void drawEffvsHitDerivedFluxSlitsCombined(TString label, TString figfile, 
					  int select_slits) {
  myCanvas->Clear(); 

  TGraphErrors *g1 = drawEffvsHitDerivedFluxSlits(label, figfile, select_slits, false); 
  
  TGraphErrors *g2 = drawEffvsHitDerivedFluxSlitsByDC(label, figfile, select_slits, false); 

  TMultiGraph *mg = new TMultiGraph();
  
  mg->Add(g1); 
  mg->Add(g2); 
    
  mg->Draw("AEP");
  mg->SetTitle(Form("slits of #pm %d", select_slits));
  mg->GetXaxis()->SetTitle("Hit-derived flux [MHz/cm^{2}]");
  mg->GetYaxis()->SetTitle("Efficiency [%]");
  mg->SetMinimum(40.);

  TLegend* leg = new TLegend(0.2, 0.2, .4, .4);
  leg->AddEntry(g1, "whole ROC"); 
  leg->AddEntry(g2, "double column"); 
  leg->SetFillColor(0);
  leg->Draw();

  myCanvas->Update(); 
  myCanvas->Print(figfile);

}

void drawEffvsHitDerivedFluxSlitsCombinedDC(TString label, TString figfile, 
					    int select_slits, int first_dc=1, 
					    int last_dc=5) {
  myCanvas->Clear(); 
  TLegend* leg = new TLegend(0.1, 0.1, .2, .2);

  TMultiGraph *mg = new TMultiGraph();
  int color = 0; 
  for (int dc = first_dc; dc <= last_dc; dc++ ) {
    color += 2; 
    TGraphErrors *g = drawEffvsHitDerivedFluxSlitsByDC(label, figfile, select_slits, false, dc, color);

    // TProfile *g = drawEffvsHitDerivedFluxSlitsByDCProfile(label, figfile, select_slits, false, dc, color);

    
    mg->Add(g);
    leg->AddEntry(g, Form("DC = %d", dc));
  } 
    
  mg->Draw("AEP");
  mg->SetTitle(Form("slits of #pm %d", select_slits));
  mg->GetXaxis()->SetTitle("Hit-derived flux [MHz/cm^{2}]");
  mg->GetXaxis()->SetLimits(0, 1000);
  mg->GetYaxis()->SetTitle("Efficiency [%]");
  mg->SetMinimum(40.);
  mg->SetMaximum(100.);

  leg->SetFillColor(0);
  leg->Draw();

  myCanvas->Update(); 
  myCanvas->Print(figfile);
  myCanvas->Write(); 

}



void draw(TString label) {
 
  TString figpath = "/afs/cern.ch/user/x/xshi/www/pxl/fig"; 
  TString figfile = Form("%s/efficiency_%s.pdf", figpath.Data(), label.Data()); 
  TString rootfile = Form("%s/efficiency_%s.root", figpath.Data(), label.Data()); 

  if ( label != "v4.7") {
    cerr << "Use label: v4.7" << label << endl; 
    return ; 
  }

  initialize("v4");

  TFile *f = new TFile(rootfile, "RECREATE");

  myCanvas = new TCanvas("c", "c", 600, 600);
  set_root_style(); 
  myCanvas->UseCurrentStyle(); 
  myCanvas->Print(Form("%s[", figfile.Data()));

  // drawEffvsHitDerivedFluxSlitsCombinedDC(label, figfile, 20);
  TH1D *g = drawEffvsHitDerivedFluxSlitsByDCProfile(label, figfile, 20, false, 5, 2);
  g->Draw(); 
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
  draw(label); 

  gSystem->Exit(0);

  return 0 ;
}

#endif

