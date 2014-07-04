#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TText.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TMath.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TStyle.h>
#include <TColor.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TObjArray.h>
#include <TMultiGraph.h>

#include "drawWithStyle.C"

#include "smoothing/interpolate.h"
#include "smoothing/SmoothingUtils.C"
#include "smoothing/LimitSmoothing.C"


Bool_t equal(double number, double reference, double tolerance = 1e-3, double epsilon = 1e-6) {
  if ((TMath::Abs(number) < epsilon) && (TMath::Abs(reference) < epsilon))
    return kTRUE;
  return TMath::Abs(number - reference) <= tolerance * TMath::Abs(number + reference);
}

using namespace std;

//_____________________________________________________________________________
class Segment {
public:
  vector<double>    x, y;
  bool              used;

  Segment(double x1, double y1, double x2, double y2)
    : x(2), y(2), used(false)
  {
    x[0]  = x1;   x[1]  = x2;
    y[0]  = y1;   y[1]  = y2;
  }
  int joins(int edge, const Segment& other) const 
  {
    for (unsigned int index = 0; index < other.x.size(); ++index) {
      if (equal(x[edge], other.x[index]) && equal(y[edge], other.y[index]))
        return index;
    } // end loop over other's endpoints
    return -1;
  }
};
//TMultiGraph* drawOutline(TH2* histogram, double threshold, Pens& pens, int index, double cornerSize = 0.25)
TMultiGraph* drawOutline(TH2* histogram, double threshold, double cornerSize = 0.25)
{
  const TAxis*                    xAxis       = histogram->GetXaxis();
  const TAxis*                    yAxis       = histogram->GetYaxis();
  const int                       numX        = xAxis->GetNbins();
  const int                       numY        = yAxis->GetNbins();

  //-- Register external edges ------------------------------------------------
  vector<Segment>                 edges;
  for (int xBin = 1; xBin <= numX; ++xBin) {
    for (int yBin = 1; yBin <= numY; ++yBin) {
      double                      content     = histogram->GetBinContent(xBin, yBin);
      if (content < threshold)    continue;
      const bool                  leftEdge    = (xBin == 1    || histogram->GetBinContent(xBin-1, yBin  ) < threshold);
      const bool                  rightEdge   = (xBin == numX || histogram->GetBinContent(xBin+1, yBin  ) < threshold);
      const bool                  bottomEdge  = (yBin == 1    || histogram->GetBinContent(xBin  , yBin-1) < threshold);
      const bool                  topEdge     = (yBin == numY || histogram->GetBinContent(xBin  , yBin+1) < threshold);
      const double                xWidth      = xAxis->GetBinWidth  (xBin);
      const double                yWidth      = yAxis->GetBinWidth  (yBin);
      double                      leftX       = xAxis->GetBinLowEdge(xBin)  , rightX      = xAxis->GetBinUpEdge (xBin);
      double                      leftYBot    = yAxis->GetBinLowEdge(yBin)  , rightYBot   = yAxis->GetBinLowEdge(yBin);
      double                      leftYTop    = yAxis->GetBinUpEdge (yBin)  , rightYTop   = yAxis->GetBinUpEdge (yBin);
      double                      bottomY     = yAxis->GetBinLowEdge(yBin)  , topY        = yAxis->GetBinUpEdge (yBin);
      double                      bottomXLft  = xAxis->GetBinLowEdge(xBin)  , topXLft     = xAxis->GetBinLowEdge(xBin);
      double                      bottomXRgt  = xAxis->GetBinUpEdge (xBin)  , topXRgt     = xAxis->GetBinUpEdge (xBin);

      //.. Register corners ...................................................
      if (bottomEdge && leftEdge  ) {
        bottomXLft               += cornerSize * xWidth;
        leftYBot                 += cornerSize * yWidth;
        edges.push_back(Segment(bottomXLft, bottomY, leftX, leftYBot));
      }
      if (leftEdge   && topEdge   ) {
        leftYTop                 -= cornerSize * yWidth;
        topXLft                  += cornerSize * xWidth;
        edges.push_back(Segment(leftX, leftYTop, topXLft, topY));
      }
      if (topEdge    && rightEdge ) {
        topXRgt                  -= cornerSize * xWidth;
        rightYTop                -= cornerSize * yWidth;
        edges.push_back(Segment(topXRgt, topY, rightX, rightYTop));
      }
      if (rightEdge  && bottomEdge) {
        rightYBot                += cornerSize * yWidth;
        bottomXRgt               -= cornerSize * xWidth;
        edges.push_back(Segment(rightX, rightYBot, bottomXRgt, bottomY));
      }

      //.. Register edges .....................................................
      if (leftEdge  )             edges.push_back(Segment(leftX     , leftYBot  , leftX     , leftYTop  ));
      if (topEdge   )             edges.push_back(Segment(topXLft   , topY      , topXRgt   , topY      ));
      if (rightEdge )             edges.push_back(Segment(rightX    , rightYBot , rightX    , rightYTop ));
      if (bottomEdge)             edges.push_back(Segment(bottomXLft, bottomY   , bottomXRgt, bottomY   ));
    } // end loop over x-bins
  } // end loop over x-bins

  //-- Link connected regions -------------------------------------------------
  TMultiGraph*                    outline     = new TMultiGraph;
  for (unsigned int iEdge = 0; iEdge < edges.size(); ++iEdge) {
    const Segment&                seed        = edges[iEdge];
    if (seed.used)                continue;

    const Segment*                link        = &seed;
    int                           attachTo    = 1;
    TGraph*                       region      = new TGraph;
    //pens.apply(region, index);
    region->SetPoint(region->GetN(), seed.x[0], seed.y[0]);
    region->SetPoint(region->GetN(), seed.x[1], seed.y[1]);

    ////cout << "_______________________________________________________________________________" << endl;
    ////cout << TString::Format("(%8.4g,%8.4g)--(%8.4g,%8.4g)", link->x[0], link->y[0], link->x[1], link->y[1]) << endl;

    for (bool canGrow = true; canGrow;) {
      canGrow                     = false;
      for (unsigned int jEdge = iEdge+1; jEdge < edges.size(); ++jEdge) {
        Segment&                  candidate   = edges[jEdge];
        if (candidate.used)       continue;
        const int                 joinedTo    = link->joins(attachTo, candidate);
        if (joinedTo >= 0) {
          region->SetPoint(region->GetN(), candidate.x[!joinedTo], candidate.y[!joinedTo]);
          candidate.used          = true;
          link                    = &candidate;
          attachTo                = !joinedTo;
          canGrow                 = true;
          ////cout << TString::Format("  <=>  %1s(%8.4g,%8.4g)--(%8.4g,%8.4g)%1s", joinedTo==0 ? "*" : "", link->x[0], link->y[0], link->x[1], link->y[1], joinedTo==1 ? "*" : "") << endl;
          break;
        }
        ////cout << TString::Format("   ?    (%8.4g,%8.4g)--(%8.4g,%8.4g) ", link->x[0], link->y[0], link->x[1], link->y[1]) << endl;
      } // end loop over candidate edges
    } // end loop over edges to link

    assert(link != &seed);
    assert(link->joins(attachTo, seed) == 0);
    region->SetPoint(region->GetN(), link->x[attachTo], link->y[attachTo]);
    outline->Add(region);
  } // end loop over unused edges 

  return outline;
}

TCanvas *pippo(TFile* file, const TString &model, const TString &scenarioX, const TString& polschema, bool pol, const TString &type, const TString &limitType){

  TH2F* hExclusion    = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_excludedPoints");
  TH2F* hExclusion_p1 = 0;
  TH2F* hExclusion_m1 = 0;
  
  if(!pol){
    hExclusion_p1 = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_p1s_excludedPoints");
    hExclusion_m1 = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_m1s_excludedPoints");
  }
  else{
    hExclusion_p1 = (TH2F*)file->Get(model + "_" + scenarioX + "_expected_R_excludedPoints");
    hExclusion_m1 = (TH2F*)file->Get(model + "_" + scenarioX + "_expected_L_excludedPoints");
  }
  
  assert(hExclusion    !=0);
  assert(hExclusion_p1 !=0);
  assert(hExclusion_m1 !=0);

  TCanvas* c2 = new TCanvas();
  TH2F* hlimit_exp = 0;
  
  if(type == "x"){ // x stays for xsection
    hlimit_exp     = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_xsection_UL");
    cout << "Min,max: " << hlimit_exp->GetMinimum(0.00001) << "," << hlimit_exp->GetMaximum() << endl;
    hlimit_exp->SetMaximum(1e2);  
    hlimit_exp->SetMinimum(2e-3);  
    c2->SetLogz(1);
    embellish(hlimit_exp, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma");
  }
  
  else if(type == "s"){ // s stays for strength
    hlimit_exp     = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_strength_UL");
    hlimit_exp->SetMaximum(2.5);
    embellish(hlimit_exp, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma/#sigma_{SUSY}");
  }
  else{
    cout << type << ": unknown type. *** Abort ***"<<endl;
    return new TCanvas();
  }
  
  TH2F *hlimit_exp_ = (TH2F*)hlimit_exp->Clone();

  hlimit_exp_->Draw("colz0text");

  TMultiGraph* outline     = drawOutline(hExclusion   , 0.999, 0.1);
  TMultiGraph* outline_p1  = drawOutline(hExclusion_p1, 0.999, 0.1);
  TMultiGraph* outline_m1  = drawOutline(hExclusion_m1, 0.999, 0.1);   
  assert(outline    !=0);
  assert(outline_p1 !=0);
  assert(outline_m1 !=0);
 
  outline    ->Draw("l");
  outline_p1 ->Draw("l");
  outline_m1 ->Draw("l");
  
  // ----------- Draw legend ----------------
  
  TLegend *leg = new TLegend(0.18,0.77,0.36,0.90);
  leg->SetTextSize(0.025);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  if(pol){
    TString name = polschema;
    if(polschema == "") name = "pol";
    c2->SetName(TString(hlimit_exp->GetTitle())+"_"+name);

    setStyle(outline   , kViolet-3, 1);
    setStyle(outline_m1, kRed,      8);
    setStyle(outline_p1, kSpring,   7);
   
    leg->AddEntry(outline->GetListOfGraphs()->First(),"Unpolarized","l");

    if(polschema == ""){
      leg->AddEntry(outline_p1->GetListOfGraphs()->First(),"Right","l");
      leg->AddEntry(outline_m1->GetListOfGraphs()->First(),"Left","l");
    }
    else{
      TObjArray*  tarray= polschema.Tokenize("_");
      TString down(tarray->At(0)->GetName());
      TString up(tarray->At(1)->GetName());
      if(down == "LL") leg->AddEntry(outline_m1->GetListOfGraphs()->First(),"#tilde{#chi}^{#pm}_{left} W_{left}","l");
      if(down == "RL") leg->AddEntry(outline_m1->GetListOfGraphs()->First(),"#tilde{#chi}^{#pm}_{right} W_{left}","l");
      if(up   == "LR") leg->AddEntry(outline_p1->GetListOfGraphs()->First(),"#tilde{#chi}^{#pm}_{left} W_{right}","l");
      if(up   == "RR") leg->AddEntry(outline_p1->GetListOfGraphs()->First(),"#tilde{#chi}^{#pm}_{right} W_{right}","l");
    }

  }
  else{
    c2->SetName(hlimit_exp->GetTitle());
  
    int col = 0;
    if(limitType == "expected") col = kRed;
    if(limitType == "observed") col = kBlack;

    setStyle(outline   , col, 1);
    setStyle(outline_m1, col, 7);
    setStyle(outline_p1, col, 7);

    leg->AddEntry(outline_p1->GetListOfGraphs()->First()," ","l");
    if(limitType == "expected") leg->AddEntry(outline->GetListOfGraphs()->First(),"Expected, #pm 1 #sigma_{experiment}","l");
    if(limitType == "observed") leg->AddEntry(outline->GetListOfGraphs()->First(),"Observed, #pm 1 #sigma_{theory}","l");
    leg->AddEntry(outline_m1->GetListOfGraphs()->First()," ","l");
  }

  hlimit_exp->SetTitle("");





  // ----------------------------------------------------------------------------------------------------------------------
  // Add expected limits on top of observed

  if(limitType == "observed"){
    TH2F* hExpExclusion    = (TH2F*)file->Get(model + "_" + scenarioX + "_expected_excludedPoints");
    TH2F* hExpExclusion_p1 = (TH2F*)file->Get(model + "_" + scenarioX + "_expected_p1s_excludedPoints");
    TH2F* hExpExclusion_m1 = (TH2F*)file->Get(model + "_" + scenarioX + "_expected_m1s_excludedPoints");
  
    assert(hExpExclusion    !=0);
    assert(hExpExclusion_p1 !=0);
    assert(hExpExclusion_m1 !=0);

    TMultiGraph* exp_outline     = drawOutline(hExpExclusion   , 0.999, 0.1);
    TMultiGraph* exp_outline_p1  = drawOutline(hExpExclusion_p1, 0.999, 0.1);
    TMultiGraph* exp_outline_m1  = drawOutline(hExpExclusion_m1, 0.999, 0.1);   
    assert(exp_outline    !=0);
    assert(exp_outline_p1 !=0);
    assert(exp_outline_m1 !=0);
 
    exp_outline    ->Draw("l");
    exp_outline_p1 ->Draw("l");
    exp_outline_m1 ->Draw("l");

    setStyle(exp_outline   , kRed, 1);
    setStyle(exp_outline_m1, kRed, 7);
    setStyle(exp_outline_p1, kRed, 7);

  
    leg->AddEntry(exp_outline_p1->GetListOfGraphs()->First()," ","l");
    leg->AddEntry(exp_outline->GetListOfGraphs()->First(),"Expected, #pm 1 #sigma_{experiment}","l");
    leg->AddEntry(exp_outline_m1->GetListOfGraphs()->First()," ","l");
  }  
  
  // ----------------------------------------------------------------------------------------------------------------------




  leg->Draw("same");
  
  // ---------------------------------------
  return c2;
}

TList *getOutline(TH2F *histo){
  TGraph2D *tg2d = new TGraph2D(histo); //Crea TGraph2D da TH2
  Double_t contours[1]; // Crea array con i valori dei contour che vuoi
  contours[0] = 1.0;
  tg2d->GetHistogram()->SetContour(1,contours);  //SetContour(numero contour,array contour)
  tg2d->Draw("cont list"); //Dummy plotting, serve solo per creare la lista di contour
  TList *contLevel = tg2d->GetContourList(1.); // Prendi quello che ti interessa
  return contLevel; 
}



TCanvas *getExclusionPlot(TFile* file, const TString &model, const TString &scenarioX, const TString& polschema, bool pol, const TString &type, const TString &limitType){

  TCanvas* c2 = new TCanvas();
  TH2F* hlimit_exp     = 0;

  if(type == "x"){ // x stays for xsection

    // ----- Getting the background for the plot, smoothing it out and embellishing it -----
    TH2F *hlimit_exp_     = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_xsection_UL");
    cout << "Min,max: " << hlimit_exp_->GetMinimum(0.00001) << "," << hlimit_exp_->GetMaximum() << endl;
    hlimit_exp = (TH2F*)hlimit_exp_->Clone();
    hlimit_exp = interpolate(hlimit_exp_,"SW"); hlimit_exp = interpolate(hlimit_exp,"SW"); 
    hlimit_exp->SetMaximum(1e2);  
    hlimit_exp->SetMinimum(2e-3);  
    c2->SetLogz(1);
    embellish(hlimit_exp, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma");
    // -------------------------------------------------------------------------------------

  }
  
  else if(type == "s"){ // s stays for strength

    // ----- Getting the background for the plot, smoothing it out and embellishing it -----
    TH2F *hlimit_exp_     = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_strength_UL");
    cout << "Min,max: " << hlimit_exp_->GetMinimum() << "," << hlimit_exp_->GetMaximum() << endl;
    hlimit_exp = (TH2F*)hlimit_exp_->Clone();
    hlimit_exp = interpolate(hlimit_exp_,"SW"); hlimit_exp = interpolate(hlimit_exp,"SW"); 
    hlimit_exp->SetMaximum(2.5);
    embellish(hlimit_exp, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma/#sigma_{SUSY}");
    // -------------------------------------------------------------------------------------

  }
  else{
    cout << type << ": unknown type. *** Abort ***"<<endl;
    return new TCanvas();
  }
 
  hlimit_exp->Draw("colz0");

















  return c2;

}


void makePlots_smoothing(TString model = "T2tt", TString scenarioX = "", TString polschema = "", bool pol = false, TString limitType = "expected"){
  TStyle* myStyle = setTDRStyle();
  paletteColdToHot(myStyle,"TChiWX");
    
  TFile *file     = new TFile(model+"_"+scenarioX+"_"+polschema+"_sigma_UL_bestexpected.root");

  // ------------------ Best region used in the limit ------------------
  TCanvas* c1 = new TCanvas();
  c1->SetName(model + "_" + scenarioX + "_bestRegion");
  TH2F* hbestRegion     = (TH2F*)file->Get(model + "_" + scenarioX +  "_bestRegion");
  embellish(hbestRegion, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]", "Region");
  hbestRegion->SetTitle("");
  hbestRegion->Draw("col text");
 
  // ----------- Draw text ----------------
  
  TLatex l;
  l.SetTextAlign(12);
  l.SetTextSize(0.035);
  l.SetTextFont(132);
  
  l.SetNDC();
  
  l.DrawLatex(0.20,0.972,"CMS Preliminary");
  l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 19.03 fb^{-1}");
  
  TString s_top;
  char *s_LSP = new char[100];
  
  if(model == "T2tt"){
    s_top = "pp#rightarrow #tilde{t} #tilde{t}*; #tilde{t}#rightarrow t + #tilde{#chi}^{0}";
    l.DrawLatex(0.57,0.965,s_top);
  }
  if(model == "T2bw"){
    s_top = "pp#rightarrow #tilde{t} #tilde{t}*; #tilde{t}#rightarrow b + #tilde{#chi}^{+}; #tilde{#chi}^{+} #rightarrow W^{+} + #tilde{#chi}^{0}";
    TString scenario = "";
    if(scenarioX == "0p25") scenario = "0.25";
    if(scenarioX == "0p50") scenario = "0.50";
    if(scenarioX == "0p75") scenario = "0.75";
    
    sprintf(s_LSP,"x = %s",scenario.Data());
    l.DrawLatex(0.63477,0.935,s_LSP);
    l.DrawLatex(0.50,0.978,s_top);
  }

  c1->SaveAs(".pdf");
  //c1->SaveAs(".root");

  // -------------------------------------------------------------------
  

  // ------------------ Excluded region --------------------------------

  
  // ------ Strength -----
  TCanvas *c2 = getExclusionPlot(file, model, scenarioX, polschema, pol, "s", limitType);
  c2->cd();

  // ----------- Draw text ----------------
  
  l.DrawLatex(0.20,0.972,"CMS Preliminary");
  l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 19.03 fb^{-1}");
  
  if(model == "T2tt")
      l.DrawLatex(0.57,0.965,s_top);
  
  if(model == "T2bw"){
    l.DrawLatex(0.63477,0.935,s_LSP);
    l.DrawLatex(0.50,0.978,s_top);
  }

  c2->SaveAs(".pdf");
  c2->SaveAs(".root");

  // Cross section UL

  TCanvas *c3 = getExclusionPlot(file, model, scenarioX, polschema, pol, "x", limitType);
  c3->cd();

  // ----------- Draw text ----------------
  
  l.DrawLatex(0.20,0.972,"CMS Preliminary");
  l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 19.03 fb^{-1}");
  
  if(model == "T2tt")
      l.DrawLatex(0.57,0.965,s_top);
  
  if(model == "T2bw"){
    l.DrawLatex(0.63477,0.935,s_LSP);
    l.DrawLatex(0.50,0.978,s_top);
  }

  c3->SaveAs(".pdf");
  c3->SaveAs(".root");
  // --------------------------------------


  return;
}





