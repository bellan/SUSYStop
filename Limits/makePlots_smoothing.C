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
#include "./TGraphSmooth.h"

#include <ctime>
#include <iostream>

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

TList *getOutline(TH2F *histo){

  TGraph2D *tg2d = new TGraph2D(histo); //Crea TGraph2D da TH2

  Double_t contours[1]; // Crea array con i valori dei contour che vuoi

  tg2d->GetHistogram()->SetContour(1,contours);  //SetContour(numero contour,array contour)

  TList *contLevel = tg2d->GetContourList(1.); // Prendi quello che ti interessa

  return contLevel; 
}


TList *drawContours(TFile *file, const TString &histoName, int col, int sty, int size){
  TH2F* hExclusion_    = (TH2F*)file->Get(histoName);
  //TH2F* hExclusion    = rebin(hExclusion_   ,"SW"); hExclusion    = rebin(hExclusion   ,"SW");
  TH2F* hExclusion    = (TH2F*)hExclusion_->Clone();
  TList *outline    = getOutline(hExclusion);
  assert(outline    !=0);

  TGraphSmooth *gs = new TGraphSmooth("normal");
  TIter next(outline);
  TObject *contour = 0;
  int ncontours = 0;
  while ((contour = next())){ // loop over the possible disjoint regions
    ++ncontours;
    setStyle((TGraph*)(contour)   , col, sty, size);

    TGraph *fGin = (TGraph*)contour;
    int fNin = fGin->GetN();
    Double_t *xinA = new Double_t[fNin];
    Double_t *yinA = new Double_t[fNin];
    Double_t *xinB = new Double_t[fNin];
    Double_t *yinB = new Double_t[fNin];

    int typeA = 0;
    int typeB = 0;
    for (int i=0;i<fNin;i++) {
      if (i != 0)
	if(fGin->GetX()[i-1] < fGin->GetX()[i]){
	  xinA[i] = fGin->GetX()[i];
	  yinA[i] = fGin->GetY()[i];
	  ++typeA;
	}
	else{
	  ++typeB;
	  xinB[i] = fGin->GetX()[i];
	  yinB[i] = fGin->GetY()[i];
	}
      else{
	++typeA;
	xinA[i] = fGin->GetX()[i];
	yinA[i] = fGin->GetY()[i];
      }
    }
    TGraph *graphA = new TGraph(typeA,xinA,yinA);
    TGraph *graphB = 0;
    if (typeB != 0) 
      graphB = new TGraph(typeB,xinB,yinB);



    //contour = gs->SmoothKern((TGraph*)contour,"normal",5.0);
    TGraph* smoothed = (TGraph*)gs->SmoothLowess((TGraph*)contour,"",0.2);
    //TGraph* smoothed = (TGraph*)gs->SmoothSuper((TGraph*)contour,"",3);
    //TGraph* smoothed = (TGraph*)contour;  

    setStyle(smoothed   , col, sty, size);
    smoothed->DrawClone("LX"); 
    
    for(unsigned i=0; i < ((TGraph*)contour)->GetN(); ++i){
      double x,y;
      ((TGraph*)contour)->GetPoint(i,x,y);
      cout << i << " x: " << x << " y: " << y << endl;
    }

  }
  cout << "N contours for " << histoName << " is " << ncontours << endl;
  return outline;
}


TCanvas *getExclusionPlot(TFile* file, const TString &model, const TString &scenarioX, const TString& polschema, bool pol, const TString &type, const TString &limitType){  

  TCanvas* c2 = new TCanvas();
  TH2F* hlimit_exp     = 0;
  TH2F * h_fake = 0; // just for binning

  //TH2F *hlimit_exp_    = 0;
  if(type == "x"){ // x stays for xsection

    // ----- Getting the background for the plot, smoothing it out and embellishing it -----
    TH2F *hlimit_exp_     = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_xsection_UL");
    cout << "Min,max: " << hlimit_exp_->GetMinimum(0.00001) << "," << hlimit_exp_->GetMaximum() << endl;
    hlimit_exp = (TH2F*)hlimit_exp_->Clone();
    //hlimit_exp = rebin(hlimit_exp,"SW"); hlimit_exp = rebin(hlimit_exp,"SW"); hlimit_exp = rebin(hlimit_exp,"SW"); hlimit_exp = rebin(hlimit_exp,"SW"); hlimit_exp = rebin(hlimit_exp,"SW"); 
    hlimit_exp->SetMaximum(1e2);  
    hlimit_exp->SetMinimum(2e-3);  
    c2->SetLogz(1);
    embellish(hlimit_exp, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma");
    // -------------------------------------------------------------------------------------
    h_fake = new TH2F("h_fake", "h_fake", 
			    hlimit_exp->GetXaxis()->GetNbins()*2, hlimit_exp->GetXaxis()->GetXmin(), hlimit_exp->GetXaxis()->GetXmax(), 
			    hlimit_exp->GetYaxis()->GetNbins()*2, hlimit_exp->GetYaxis()->GetXmin(), hlimit_exp->GetYaxis()->GetXmax()
			    );
    embellish(h_fake, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma");
  }
  
  else if(type == "s"){ // s stays for strength

    // ----- Getting the background for the plot, smoothing it out and embellishing it -----
    TH2F *hlimit_exp_     = (TH2F*)file->Get(model + "_" + scenarioX + "_"+limitType+"_strength_UL");
    cout << "Min,max: " << hlimit_exp_->GetMinimum() << "," << hlimit_exp_->GetMaximum() << endl;
    hlimit_exp = (TH2F*)hlimit_exp_->Clone();
    //hlimit_exp = interpolate(hlimit_exp_,"SW"); hlimit_exp = interpolate(hlimit_exp,"SW"); hlimit_exp = rebin(hlimit_exp,"SW"); //hlimit_exp = rebin(hlimit_exp,"SW"); 
    hlimit_exp->SetMaximum(2.5);
    embellish(hlimit_exp, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma/#sigma_{SUSY}");
    // -------------------------------------------------------------------------------------
    h_fake = new TH2F("h_fake", "h_fake", 
		      hlimit_exp->GetXaxis()->GetNbins()*2, hlimit_exp->GetXaxis()->GetXmin(), hlimit_exp->GetXaxis()->GetXmax(), 
		      hlimit_exp->GetYaxis()->GetNbins()*2, hlimit_exp->GetYaxis()->GetXmin(), hlimit_exp->GetYaxis()->GetXmax()
		      );
    embellish(h_fake, "m_{#tilde{t}} [GeV]", "m_{LSP} [GeV]","95% CL limit on #sigma/#sigma_{SUSY}"); 
  }
  else{
    cout << type << ": unknown type. *** Abort ***"<<endl;
    return new TCanvas();
  }
 
  h_fake->GetYaxis()->SetRangeUser(0,400);                                      h_fake->GetXaxis()->SetRangeUser(200, model == "T2tt" ? 900 : 800);
  h_fake->GetXaxis()->SetLabelSize(hlimit_exp->GetXaxis()->GetLabelSize());     h_fake->GetYaxis()->SetLabelSize(hlimit_exp->GetYaxis()->GetLabelSize());
  h_fake->GetXaxis()->SetLabelOffset(hlimit_exp->GetXaxis()->GetLabelOffset()); h_fake->GetYaxis()->SetLabelOffset(hlimit_exp->GetYaxis()->GetLabelOffset());
  h_fake->GetXaxis()->SetTitleSize(hlimit_exp->GetXaxis()->GetTitleSize());     h_fake->GetYaxis()->SetTitleSize(hlimit_exp->GetYaxis()->GetTitleSize());
  h_fake->GetXaxis()->SetTitleOffset(hlimit_exp->GetXaxis()->GetTitleOffset()); h_fake->GetYaxis()->SetTitleOffset(hlimit_exp->GetYaxis()->GetTitleOffset());
  
  h_fake->Draw();

  hlimit_exp->Draw("same colz0");
  //hlimit_exp->Draw("colz0");

  h_fake->Draw("same axis");
  c2->SetName(hlimit_exp->GetTitle());
  

  // Now get the normalized (to susy reference cross section) plots and draw the contours
  int col = 0;
  int sty = 1;
  if(limitType == "expected") {col = kRed; sty = 7;}
  if(limitType == "observed") col = kBlack;

  // ----------- Draw legend ----------------
  
  TLegend *leg = new TLegend(0.18,0.77,0.36,0.90);
  leg->SetTextSize(0.025);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  TString legendtag = "";
  if(limitType == "expected") legendtag = "Expected, #pm 1 #sigma_{experiment}";  
  if(limitType == "observed") legendtag = "Observed, #pm 1 #sigma_{theory}";

  leg->AddEntry((TGraph*)(drawContours(file, model + "_" + scenarioX + "_"+limitType+"_p1s_strength_UL", col, sty, 2)->First())," ","l");
  leg->AddEntry((TGraph*)(drawContours(file, model + "_" + scenarioX + "_"+limitType+"_strength_UL"    , col, sty, 4)->First()),legendtag,"l");
  leg->AddEntry((TGraph*)(drawContours(file, model + "_" + scenarioX + "_"+limitType+"_m1s_strength_UL", col, sty, 2)->First())," ","l");

  hlimit_exp->SetTitle("");

  
  // ----------------------------------------------------------------------------------------------------------------------
  // Add expected limits on top of observed
  
  if(limitType == "observed"){

    // Now get the normalized (to susy reference cross section) plots
    leg->AddEntry((TGraph2D*)(drawContours(file, model + "_" + scenarioX + "_expected_p1s_strength_UL", kRed, 7, 2)->First())," ","l");
    leg->AddEntry((TGraph2D*)(drawContours(file, model + "_" + scenarioX + "_expected_strength_UL"    , kRed, 7, 4)->First()),"Expected, #pm 1 #sigma_{experiment}","l");
    leg->AddEntry((TGraph2D*)(drawContours(file, model + "_" + scenarioX + "_expected_m1s_strength_UL", kRed, 7, 2)->First())," ","l");
  }


  leg->Draw("same");
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
  
  l.DrawLatex(0.20,0.972,"CMS");
  l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 18.9 fb^{-1}");
  
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

  
  // Cross section UL
  TCanvas *c3 = getExclusionPlot(file, model, scenarioX, polschema, pol, "x", limitType);
  c3->cd();

  // ----------- Draw text ----------------
  
  l.DrawLatex(0.20,0.972,"CMS");
  l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 18.9 fb^{-1}");
  
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

  // ------ Strength -----
  TCanvas *c2 = getExclusionPlot(file, model, scenarioX, polschema, pol, "s", limitType);
  c2->cd();

  // ----------- Draw text ----------------
  
  l.DrawLatex(0.20,0.972,"CMS");
  l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 18.9 fb^{-1}");
  
  if(model == "T2tt")
      l.DrawLatex(0.57,0.965,s_top);
  
  if(model == "T2bw"){
    l.DrawLatex(0.63477,0.935,s_LSP);
    l.DrawLatex(0.50,0.978,s_top);
  }

  c2->SaveAs(".pdf");
  c2->SaveAs(".root");


  return;
}





