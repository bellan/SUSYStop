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
//#include "./myTGraphSmooth.h"
#include <TGraphSmooth.h>

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

//class SmoothCurve{
//  SmoothCurve();
//};

std::pair<int,int> split(int fNin, const Double_t *xin, const Double_t *yin, Double_t *xoutA, Double_t *youtA, Double_t *xoutB, Double_t *youtB){
  double max = *std::max_element(xin, xin + fNin); 
  //cout << "MAX " << max <<endl;
    int typeA = 0;
    int typeB = 0;

    bool turningPoint = false;
    for (int i=0;i<fNin;++i) {
      if(xin[i] <= max && (!turningPoint || (max - xin[i]) < 5) && typeB == 0){
	xoutA[typeA] = xin[i];
	youtA[typeA] = yin[i];
	++typeA;
	//cout << "A: " << i << " x: " << xoutA[i] << " y: " << youtA[i] << endl;
      }
      else{
	if(typeB == 0){
	  xoutB[typeB] = xin[i-1];
	  youtB[typeB] = yin[i-1];
	  //cout << "B: "<< i << " " << typeB << " x: " << xoutB[typeB] << " y: " << youtB[typeB] << endl;
	  ++typeB;
	}
	xoutB[typeB] = xin[i];
	youtB[typeB] = yin[i];
	//cout << "B: "<< i << " " << typeB << " x: " << xoutB[typeB] << " y: " << youtB[typeB] << endl;
	++typeB;
      }
      if(!turningPoint && xin[i] == max) turningPoint = true;
    }
    return std::make_pair(typeA,typeB);
}


TList *drawContours(TFile *file, const TString &histoName, int col, int sty, int size){
  TH2F* hExclusion_    = (TH2F*)file->Get(histoName);
  TH2F* hExclusion    = rebin(hExclusion_   ,"SW"); //hExclusion    = rebin(hExclusion   ,"SW");
  //TH2F* hExclusion    = (TH2F*)hExclusion_->Clone();
  TList *outline    = getOutline(hExclusion);
  assert(outline    !=0);

  TIter next(outline);
  TObject *contour = 0;
  int ncontours = 0;
  cout << histoName << " ncontours: " << ncontours <<endl;
  while ((contour = next())){ // loop over the possible disjoint regions
    TGraph *fGin = (TGraph*)contour;
    int fNin = fGin->GetN();
    
    cout << "Contour: " << ncontours << " Npoints: " << fNin  << endl;
    ++ncontours;
    setStyle((TGraph*)(contour)   , col, sty, size);


    if(fNin < 10) {cout << "Rejected: " << fNin << endl; continue;}

    Double_t *xoutA = new Double_t[fNin+1];
    Double_t *youtA = new Double_t[fNin+1];
    Double_t *xoutB = new Double_t[fNin];
    Double_t *youtB = new Double_t[fNin];

    int typeA = 0;
    int typeB = 0;

    Double_t * xin = fGin->GetX();
    Double_t * yin = fGin->GetY();

    std::pair<int,int> points = split(fNin,xin,yin,xoutA,youtA,xoutB,youtB);
    typeA = points.first;
    typeB = points.second;

    //cout << "Points division " << fNin << " " << typeA << " " << typeB << endl; 

    TGraphSmooth *gs = new TGraphSmooth("normal");

    Double_t *xoutAB = new Double_t[fNin+3];
    Double_t *youtAB = new Double_t[fNin+3];

    TGraph *graphA = new TGraph(typeA,xoutA,youtA);
    //TGraph* smoothedA = graphA;
    TGraph* smoothedA = (TGraph*)gs->SmoothLowess(graphA, "",       0.2, 3, 12.5);
    setStyle(smoothedA   , kBlue, sty, size);
    //smoothedA->DrawClone("LX"); 

    xoutAB[0] = xoutA[0];
    youtAB[0] = youtA[0];

    int nAB = smoothedA->GetN()+1;
    
    Double_t *xsmoothedA = smoothedA->GetX();
    Double_t *ysmoothedA = smoothedA->GetY();

    
    for(int i=0; i < smoothedA->GetN(); ++i){
      xoutAB[i+1] = xsmoothedA[i];
      youtAB[i+1] = ysmoothedA[i];
    }


    if((xoutA[graphA->GetN()-1] != xoutAB[smoothedA->GetN()] || youtA[graphA->GetN()-1] != youtAB[smoothedA->GetN()])){
      xoutAB[smoothedA->GetN()+1] = xoutA[graphA->GetN()-1];
      youtAB[smoothedA->GetN()+1] = youtA[graphA->GetN()-1];
      ++nAB;
    }



    if (typeB != 0){ 
      TGraph *graphB = new TGraph(typeB,xoutB,youtB);
      //TGraph* smoothedB = graphB;
      TGraph* smoothedB = (TGraph*)gs->SmoothLowess(graphB, "",       0.2, 3, 12.5);
      setStyle(smoothedB   , kGreen, sty, size);
      //smoothedB->DrawClone("LX"); 
      if(graphB) delete graphB;
      
      Double_t *xsmoothedB = smoothedB->GetX();
      Double_t *ysmoothedB = smoothedB->GetY();

      for(int i=0; i < smoothedB->GetN(); ++i){
	xoutAB[i+nAB] = xsmoothedB[smoothedB->GetN()-i-1];
	youtAB[i+nAB] = ysmoothedB[smoothedB->GetN()-i-1];
      }
      nAB += smoothedB->GetN();
  }
    if(graphA) delete graphA;
    
    //for(int i=0; i<typeB; ++i){
    //  cout << i << endl;
    //  xoutA[i+typeA] = xoutB[i];
    //  youtA[i+typeA] = youtB[i];
    //}
    
    cout << "------------------------------------------------------------------" << endl;
    cout<< histoName << " contour number: " << ncontours << endl;
    int nAB_clean=0;
    Double_t *xoutAB_clean = new Double_t[nAB];
    Double_t *youtAB_clean = new Double_t[nAB];
    for(int i=0;i<nAB;++i){
      cout << "Smoothed: " << i << " " << xoutAB[i] << " " << youtAB[i] << endl;
      
      bool refine = true;

      if(refine){
	if(histoName == "T2bw_0p75_expected_m1s_strength_UL" && ncontours == 2)
	  if(i == 0) continue;
	if(histoName == "T2bw_0p50_expected_m1s_strength_UL" && ncontours == 2){
	  if(i > 129 && i < 145) {cout << "removing " << i << endl; continue;}
	  if(i == 145){
	    cout << "replacing " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[129];
	    youtAB_clean[nAB_clean] = 0;
	    ++nAB_clean;
	    continue;
	  }
	}
	if(histoName == "T2bw_0p50_observed_p1s_strength_UL" && ncontours == 2)
	  if(i == 115 || i == 116 || i == 117) continue;
	
	if(histoName == "T2bw_0p50_expected_strength_UL" && ncontours == 3)
	  if(i == 113 || i == 114 || i == 115 || i == 116 || i == 117) continue;
	
	if(histoName == "T2bw_0p50_observed_strength_UL" && ncontours ==2){
	  if(i>119 && i < 131 || i == 115 || i == 116 || i == 117) continue;
	  if(i == 131){
	    cout << "replacing " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[120];
	    youtAB_clean[nAB_clean] = 0;
	    ++nAB_clean;
	    continue;
	  }
	}
	if(histoName == "T2bw_0p50_expected_m1s_strength_UL" && ncontours==1)
	  if(i >3 && i<8) continue;
	
	if(histoName == "T2bw_0p50_observed_m1s_strength_UL" && ncontours==1){
	  if(i==0) continue;
	  if(i==1){
	    cout << "replacing " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[1];
	    youtAB_clean[nAB_clean] = (youtAB[1]+youtAB[2])/2.;
	    ++nAB_clean;
	    continue;
	  }
	}
	
	if(histoName == "T2tt__observed_p1s_strength_UL" && ncontours==1){
	  if(i==16) continue;
	  if(i<16 && i > 0){
	    cout << "Sorting " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[16-i];
	    youtAB_clean[nAB_clean] = youtAB[16-i];
	    ++nAB_clean;
	    continue;
	  }
	}
	
	if(histoName == "T2tt__observed_m1s_strength_UL" && ncontours==1){
	  if(i<22 && i > 0){
	    cout << "Sorting " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[22-i];
	    youtAB_clean[nAB_clean] = youtAB[22-i];
	    ++nAB_clean;
	    continue;
	  }
	}
	
	if(histoName == "T2tt__expected_p1s_strength_UL" && ncontours==1){
	  if(i==18) continue;
	  if(i<18 && i > 0){
	    cout << "Sorting " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[18-i];
	    youtAB_clean[nAB_clean] = youtAB[18-i];
	    ++nAB_clean;
	    continue;
	  }
	}
	if(histoName == "T2tt__expected_strength_UL" && ncontours==1){
	  if(i==16) continue;
	  if(i<15 && i > 0){
	    cout << "Sorting " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[15-i];
	    youtAB_clean[nAB_clean] = youtAB[15-i];
	    ++nAB_clean;
	    continue;
	  }
	}
	if(histoName == "T2tt__expected_m1s_strength_UL" && ncontours==1){
	  if(i==19) continue;
	  if(i<19 && i > 0){
	    cout << "Sorting " << i << endl;
	    xoutAB_clean[nAB_clean] = xoutAB[19-i];
	    youtAB_clean[nAB_clean] = youtAB[19-i];
	    ++nAB_clean;
	    continue;
	  }
	}

	if(histoName=="T2bw_0p25_observed_strength_UL" && ncontours==1){
	  if(i == 0 || i ==1 || i == 2 || i==21 || i ==22) continue;
	  if(i==3){
	    xoutAB_clean[nAB_clean] = 200;
	    youtAB_clean[nAB_clean] = youtAB[i];
	    ++nAB_clean;
	    continue;
	  }
	}

	if(histoName=="T2bw_0p25_observed_strength_UL" && ncontours==2){
	  if(i == 0){
	    xoutAB_clean[nAB_clean] = xoutAB[175];
	    youtAB_clean[nAB_clean] = youtAB[175];
	    ++nAB_clean;
	    continue;	    
	  }
	  if(i == 1 || i == 18 || (i>3 && i<11) || i ==35 || i == 37) continue;
	  if(i>0 && i <35){
	    xoutAB_clean[nAB_clean] = fGin->GetX()[i];
	    youtAB_clean[nAB_clean] = fGin->GetY()[i];
	    ++nAB_clean;
	    continue;
	  }
	}
	if(histoName=="T2bw_0p25_observed_p1s_strength_UL" && ncontours==1){
	  if(i == 0 || i == 1 || i == 2 || i == 23) continue;
	  if(i==3){
	    xoutAB_clean[nAB_clean] = 200;
	    youtAB_clean[nAB_clean] = youtAB[i];
	    ++nAB_clean;
	    continue;
	  }

	}

	if(histoName=="T2bw_0p25_observed_p1s_strength_UL" && ncontours==2)
	  if(i == 0) continue;

	if(histoName=="T2bw_0p25_observed_m1s_strength_UL" && ncontours==1){
	  if(i == 0 || i ==1 || i == 2) continue;
	  if(i==3){
	    xoutAB_clean[nAB_clean] = 200;
	    youtAB_clean[nAB_clean] = youtAB[i];
	    ++nAB_clean;
	    continue;
	  }  
	}

	if(histoName=="T2bw_0p25_expected_strength_UL" && ncontours==1){
	  if(i==0 || i==34 || i==35) continue;
	  if(i==1){
	    xoutAB_clean[nAB_clean] = 200;
	    youtAB_clean[nAB_clean] = youtAB[i];
	    ++nAB_clean;
	    continue;
	  }
	}

	if(histoName=="T2bw_0p25_expected_p1s_strength_UL" && ncontours == 7)
	  if(i >72 && i <78) continue;
	
	

	if(histoName=="T2bw_0p25_expected_m1s_strength_UL" && ncontours == 1){
	  if(i == 0){
	    xoutAB_clean[nAB_clean] = 200;
	    youtAB_clean[nAB_clean] = youtAB[i];
	    ++nAB_clean;
	    continue;
	  }
	  if(i > 18 && i < 25) continue;
	  if(i == 25) { // link to the other stub   
	    xoutAB_clean[nAB_clean] = 311.271;
	    youtAB_clean[nAB_clean] = 188.558;
	    ++nAB_clean;
	    continue;
	  }
	}
	
	if(histoName=="T2bw_0p25_expected_m1s_strength_UL" && ncontours == 4)
	  if(i > 110 && i < 117) continue;
      
	if(histoName == "T2bw_0p25_expected_m1s_strength_UL" && ncontours == 2){
	  if(i>13 && i < 22)  continue;
	  if(i == 23) { // link to the other stub  
	    xoutAB_clean[nAB_clean] = 275;
	    youtAB_clean[nAB_clean] = 36.0118;
	    ++nAB_clean;
	    continue;
	  }
	}
	if(histoName == "T2bw_0p25_expected_m1s_strength_UL" && ncontours == 6){
	  if(i == 0 || i == 1) continue;
	  if(i == 2){
	    xoutAB_clean[nAB_clean] = xoutAB[i];
	    youtAB_clean[nAB_clean] = 0;
	    ++nAB_clean;
	    continue;
	  }
	  if(i>16 && i < 24) continue;
	}





      }
      
      
      xoutAB_clean[nAB_clean] = xoutAB[i];
      youtAB_clean[nAB_clean] = youtAB[i];
      ++nAB_clean;
    }
    cout << "------------------------------------------------------------------" << endl;

    for(int i=0;i<nAB_clean;++i) cout << "Final order: " << i << " " << xoutAB_clean[i] << " " << youtAB_clean[i] << endl;
  
    TGraph* smoothedAB = new TGraph(nAB_clean,xoutAB_clean,youtAB_clean);
    //TGraph* smoothedAB = new TGraph(nAB,xoutAB,youtAB);
    setStyle(smoothedAB   , col, sty, size);
    smoothedAB->DrawClone("LX"); 

    //TGraph* smoothed = (TGraph*)gs->SmoothKern  ((TGraph*)contour, "normal", 5.0);
    //TGraph* smoothed = (TGraph*)gs->SmoothLowess((TGraph*)contour, "",       0.2);
    //TGraph* smoothed = (TGraph*)gs->SmoothSuper ((TGraph*)contour, "",       3);
    
    // uncomment to have the original curve superimposed to the smoothed one.
    //TGraph* original = (TGraph*)contour;  
    //setStyle(original   , kBlack, sty, 1);
    //for(int j=0; j < original->GetN(); ++j) cout << "Original " << j << "  " << original->GetX()[j] << " " << original->GetY()[j] << endl;
    //original->DrawClone("LX"); 
    // --------------------

    //for(unsigned i=0; i < ((TGraph*)contour)->GetN(); ++i){
    //  double x,y;
    //  ((TGraph*)contour)->GetPoint(i,x,y);
    //  cout << i << " x: " << x << " y: " << y << endl;
    //}
    
  }
  //cout << "N contours for " << histoName << " is " << ncontours << endl;
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
  
  //TLegend *leg = new TLegend(0.18,0.77,0.36,0.90);
  //leg->SetTextSize(0.025);
  TLegend *leg = new TLegend(0.18,0.73,0.36,0.86);
  leg->SetTextSize(0.0325);
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
    l.DrawLatex(0.20,0.965,s_top);
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

  //l.DrawLatex(0.20,0.972,"CMS");
  //l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 18.9 fb^{-1}");
  
  l.SetTextSize(0.04);
  l.DrawLatex(0.27,0.88,"CMS #sqrt{s} = 8 TeV   L = 18.9 fb^{-1}");
  l.SetTextSize(0.05);

  if(model == "T2tt")
    //l.DrawLatex(0.57,0.965,s_top);
    l.DrawLatex(0.33,0.965,s_top);

  if(model == "T2bw"){
    //l.DrawLatex(0.63477,0.935,s_LSP);
    //l.DrawLatex(0.50,0.978,s_top);
    l.DrawLatex(0.1,0.965,s_top);
    l.DrawLatex(0.75,0.965,s_LSP);

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





