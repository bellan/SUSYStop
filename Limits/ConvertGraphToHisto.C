#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TText.h"
//#include "TLatex.h"
//#include "THStack.h"
//#include "TDirectory.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph2D.h"
//#include "TMarker.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TStyle.h>
#include <TColor.h>
#include <TSystem.h>
#include <TLatex.h>

#include "/home/bellan/Workspace/PhenomCorrection/setTDRStyle.C"

TH1F* convert(TGraph *pGraph){
  // takes data from a graph, determines binning and fills data into histogram
  Int_t NPoints = pGraph->GetN();
  Double_t BinLimits[NPoints+1];
  // sort graph
  pGraph->Sort();
  // determine lower limit of histogram: half the distance to next point
  Double_t x0,x1,y;
  pGraph->GetPoint(0,x0,y);
  pGraph->GetPoint(1,x1,y);
  Double_t Distance = TMath::Abs(x0-x1);
  BinLimits[0] = x0 - Distance/2.;
  // now set upper limits for all the other points
  for (Int_t k = 0 ; k<NPoints-1;k++){
    pGraph->GetPoint(k,x0,y);
    pGraph->GetPoint(k+1,x1,y);
    Distance = TMath::Abs(x0-x1);
    BinLimits[k+1] = x0 + Distance/2.;}
  // for the last point set upper limit similar to first point:
  pGraph->GetPoint(NPoints-2,x0,y);
  pGraph->GetPoint(NPoints-1,x1,y);
  Distance = TMath::Abs(x0-x1);
  BinLimits[NPoints] = x1 + Distance/2.;
  // now we know the binning and can create the histogram:
  TString Name = "deltar_PLOT_eff";
  TH1F *ThisHist = new TH1F(Name,"",NPoints,BinLimits);
  ThisHist->SetOption("P");
  // now fill the histogram
  for (Int_t i = 0; i<pGraph->GetN();i++){
    Double_t x,y;
    pGraph->GetPoint(i,x,y);
    //ThisHist->SetBinContent(i+1,y);
    ThisHist->SetBinContent(i+1,100*pGraph->GetErrorY(i)/y);
    //ThisHist->SetBinError(i+1,std::min(pGraph->GetErrorY(i),1-y));
  }
  return ThisHist;
} 

void ConvertGraphToHisto(){
  string fileName="sms_referenceXSecs_8TeV.root";
  TFile *_file = TFile::Open((char *) fileName.c_str());
  TGraph *gr = (TGraph*)_file->Get("stop");
  TH1F *central = convert(gr);
  TFile out("SMS_XSECS_8TeV.root","UPDATE");
  out.cd();
  central->Write("stop_xsec_unc");
  out.Close();

}
