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

#include "drawWithStyle.C"

std::pair<TH2F*,TH2F*> getExcludedPoints(TH2F* h2dXSecUL, int sigma = 0, double referenceXSecStrength = 1., bool divideByReference = true);

using namespace std;


void addExcludedRegions(TString model = "T2tt", TString scenarioX = "", TString pol = "", bool normalize = true){

  TFile *file  = new TFile(model+"_"+scenarioX+"_"+pol+"_sigma_UL_bestexpected.root","UPDATE");
  
  TString variants[] = {"","p1s","m1s","L","R"};
  int size = sizeof(variants)/sizeof(TString);
  
  for(int i=0; i<size; ++i){

    TString input = "_expected_";
    if (variants[i] != "") input = input + variants[i] + "_";
    
    // Get the expected UL xsection
    TH2F* hlimit_exp     = (TH2F*)file->Get(model + "_" + scenarioX + input + "xsection_UL");

    // get the TGraph corresponding to the expected exclusion. 1. stands for the ref_xsec multiplier
    pair<TH2F*,TH2F*> exclusion = getExcludedPoints(hlimit_exp,  0., 1., normalize);

    file->cd();
    embellish(exclusion.first, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]");
    exclusion.first->Write(model+"_"+scenarioX+input+"excludedPoints");
    embellish(exclusion.second, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]");
    exclusion.second->SetTitle(model+"_"+scenarioX+input+"strength_UL");
    exclusion.second->Write(model+"_"+scenarioX+input+"strength_UL");
  }

  file->Close();
}


std::pair<TH2F*,TH2F*> getExcludedPoints(TH2F* h2dXSecUL, int sigma, double referenceXSecStrength, bool divideByReference) {

  string mass="stop_xsec";

  string fileName="xsections/SMS_XSECS_8TeV.root";

  // Get the reference cross section and the uncertainty on it
  TFile *_file = TFile::Open((char *) fileName.c_str());
  TH1 *hReferenceXSection     = (TH1F*) _file->Get((char *) mass.c_str());
  TH1 *hReferenceXSection_unc = (TH1F*) _file->Get((mass+"_unc").c_str());

  // convert them into graph to benefit from extrapolation between points
  TGraph *grReferenceXSection     = new TGraph(hReferenceXSection);
  TGraph *grReferenceXSection_unc = new TGraph(hReferenceXSection_unc);

  // Get the properties of the UL in input
  Int_t  nBinX = h2dXSecUL->GetXaxis()->GetNbins();
  double xLow  = h2dXSecUL->GetXaxis()->GetBinLowEdge(1), xHigh = h2dXSecUL->GetXaxis()->GetBinLowEdge(nBinX+1);

  Int_t  nBinY = h2dXSecUL->GetYaxis()->GetNbins();
  double yLow  = h2dXSecUL->GetYaxis()->GetBinLowEdge(1), yHigh = h2dXSecUL->GetYaxis()->GetBinLowEdge(nBinY+1);
  
  // Bare 2D-histo (m_LSP vs m_stop) of the reference xsection  
  TH2F * h2dReferenceXSection = new TH2F("ReferenceXSection", "ReferenceXSection", nBinX, xLow ,xHigh , nBinY , yLow , yHigh);
  TH2F * h2dExcludedPoints    = new TH2F("ExcludedPoints"   , "ExcludedPoints"   , nBinX, xLow ,xHigh , nBinY , yLow , yHigh);

  // loop over the x,y bin of the UL xsection 
  for(int i=1; i <= nBinX; i++) {
    
    // Get mstop, the uncertainty on the xsection and the xsection itself. All this quantities are mlsp independent
    double mstop           = h2dXSecUL->GetXaxis()->GetBinCenter(i);
    double xsecunc         = grReferenceXSection_unc->Eval(mstop)/100.;
    double reference_xsec  = referenceXSecStrength*(1+sigma*xsecunc)*grReferenceXSection->Eval(mstop);  
    
    for(int j=1; j <= nBinY; j++){
      double mlsp          = h2dXSecUL->GetYaxis()->GetBinCenter(j);
      // Fill the 2D reference xsec (use the same code for +/-sigma_theo).
      if(h2dXSecUL->GetBinContent(i,j) != 0) h2dReferenceXSection->Fill(mstop, mlsp, reference_xsec);
    }
  }
  
  // Divide, or not, the UL xsec by the reference. The idea is then to check against unity if the point is excluded or not.
  TH2F* h2dXSecUL_OverReference = (TH2F*) h2dXSecUL->Clone("XSecUL_OverReference");
  if(divideByReference) h2dXSecUL_OverReference->Divide(h2dReferenceXSection);
  
  // Loop over x,y points of the plane
  for(int i=1; i <= nBinX; i++) {
    double mstop = h2dXSecUL_OverReference->GetXaxis()->GetBinCenter(i); 
    std::vector<double> excludedPoints;
    for(int j=1; j <= nBinY; j++) {

      double mlsp   = h2dXSecUL_OverReference->GetYaxis()->GetBinCenter(j);
      double xsecUL = h2dXSecUL_OverReference->GetBinContent(i,j);

      double upperXsec_reference = 1.; // if sigma_UL is normalized to the reference crossection
      // If the UL_xsec is not normalized to the reference, then use the xsec_ref as threasold
      if(!divideByReference){
	// use the original TH1F because it has more bins
	double xsecunc      = grReferenceXSection_unc->Eval(mstop)/100.;
	upperXsec_reference = referenceXSecStrength*(1+sigma*xsecunc)*grReferenceXSection->Eval(mstop); 
      }

      // Check if the (m_stop,m_LSP) is excluded
      if(xsecUL > 0. && xsecUL <= upperXsec_reference) {
	//cout << "Exclude m_LSP = " << mlsp << " because " << xsecUL << " <= " << upperXsec_reference << endl;
	excludedPoints.push_back(mlsp);
	h2dExcludedPoints->SetBinContent(i,j,1);
      }
      else h2dExcludedPoints->SetBinContent(i,j,0);
    }
  }

  return std::make_pair(h2dExcludedPoints,h2dXSecUL_OverReference);
}

