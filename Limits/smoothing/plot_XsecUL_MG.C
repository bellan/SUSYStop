#include "SMS_utils.C"
#include "setTDRStyle.h"
#include <stdlib.h> 
#include "interpolate.h"
#include "SmoothingUtils.C"
#include "LimitSmoothing.C"

void plot_XsecUL_MG(TString path2file, int whichHT, std::string whichSample, bool saveOutputFile = false, TString calcWhat, bool toys,int thUnc = 0) {

  // calcWhat options
  // Obs , Exp , Exp_plus1sigma , Exp_minus1sigma , Obs_thUp , Obs_thDown

  setTDRStyle();
  gStyle->SetPalette(1);

  bool fillEveryBin = false;

  TH2F *h_xSec_limit= new TH2F("h_xSec_limit","xSec_limit",81,-12.5,2012.5,81,-12.5,2012.5);


  std::ifstream file(path2file);
  std::string line;

  double num[6];
  double xsec = 0.1;

  while (std::getline(file, line)) {
    
    istringstream ss(line);
    double d;
    
    int i=0;
    while (ss >> d) { 
      num[i]=d;
      i++;
    }

    //    std::cout << num[0] << "\t" << num[1] << "\t" << num[2] << std::endl;
    
    // Fill in the 2D Histogram
    if (!toys) {
      if ( (calcWhat=="Exp") && (num[2]>0.) && (num[2]<998.)) { h_xSec_limit->SetBinContent((num[0]/25),(num[1]/25),num[2]*xsec); }
      if ( (calcWhat=="ExpUp") && (num[3]>0.) && (num[3]<998.)) { h_xSec_limit->SetBinContent((num[0]/25),(num[1]/25),num[3]*xsec); }
      if ( (calcWhat=="ExpDown") && (num[4]>0.) && (num[4]<998.)) { h_xSec_limit->SetBinContent((num[0]/25),(num[1]/25),num[4]*xsec); }
      if ( (calcWhat=="Obs") && (num[5]>0.) && (num[5]<998.)) { h_xSec_limit->SetBinContent((num[0]/25),(num[1]/25),num[5]*xsec); }
    }

    //    else if ((num[2]>0.) && (num[2]<998.) && (num[1]>0)) { h_xSec_limit->SetBinContent((num[0]/25)+1,(num[1]/25)+1,num[2]*xsec); }
    else if ((num[2]>0.) && (num[2]<998.)) { h_xSec_limit->SetBinContent((num[0]/25)+1,(num[1]/25)+1,num[2]*xsec); }
  }
  
  TH2F *h_xSec_limit_ = (TH2F*)h_xSec_limit->Clone("h_xSec_limit_");


  if (whichSample == "T1tttt") { h_xSec_limit_ = interpolate(h_xSec_limit,"SW"); h_xSec_limit = interpolate(h_xSec_limit,"SW"); } //h_xSec_limit = rebin(h_xSec_limit,"SW");
  //    h_xSec_limit = rebin(h_xSec_limit,"SW"); h_xSec_limit = rebin(h_xSec_limit,"SW"); h_xSec_limit = rebin(h_xSec_limit,"SW"); h_xSec_limit = rebin(h_xSec_limit,"SW"); h_xSec_limit = doSmooth(h_xSec_limit,true); }
  if (whichSample == "T1t1t")  { h_xSec_limit_ = interpolate(h_xSec_limit,"EW"); h_xSec_limit = interpolate(h_xSec_limit,"EW"); h_xSec_limit = rebin(h_xSec_limit,"EW");
    h_xSec_limit = rebin(h_xSec_limit,"EW"); h_xSec_limit = rebin(h_xSec_limit,"EW"); h_xSec_limit = rebin(h_xSec_limit,"EW"); h_xSec_limit = rebin(h_xSec_limit,"EW"); }
  if (whichSample == "T5tttt") { h_xSec_limit_ = interpolate(h_xSec_limit,"SW"); h_xSec_limit = interpolate(h_xSec_limit,"SW"); h_xSec_limit = rebin(h_xSec_limit,"SW"); 
    h_xSec_limit = rebin(h_xSec_limit,"SW"); h_xSec_limit = rebin(h_xSec_limit,"SW"); h_xSec_limit = rebin(h_xSec_limit,"SW"); h_xSec_limit = rebin(h_xSec_limit,"SW"); }


  if (fillEveryBin) {

    for (UInt_t iy=1; iy<h_xSec_limit->GetNbinsY()+1; iy++) {
      for (UInt_t ix=1; ix<h_xSec_limit->GetNbinsX()+1; ix++) {

	double Mx = h_xSec_limit->GetBinCenter(ix);
	double My = h_xSec_limit->GetBinCenter(iy);
	
	if ((My<=(Mx-225.)) && (h_xSec_limit->GetBinContent(ix,iy)==0) && (Mx>=400.) && (Mx<=1400.)) {

	  double previousBinX = h_xSec_limit->GetBinContent(ix-1,iy);
	  double nextBinX = h_xSec_limit->GetBinContent(ix+1,iy);
	  double previousBinY = h_xSec_limit->GetBinContent(ix,iy-1);
	  double nextBinY = h_xSec_limit->GetBinContent(ix,iy+1);
	  double S4_UL = h_xSec_limit->GetBinContent(ix-1,iy+1);
	  double S4_UR = h_xSec_limit->GetBinContent(ix+1,iy+1);
	  double S4_DL = h_xSec_limit->GetBinContent(ix-1,iy-1);
	  double S4_DR = h_xSec_limit->GetBinContent(ix+1,iy-1);

	  double mean = 0.;
	  if ((previousBinX>0.) && (nextBinX>0.)) { mean = (previousBinX + nextBinX)/2.; }
	  else if ((previousBinY>0.) && (nextBinY>0.)) { mean = (previousBinY + nextBinY)/2.; }
	  else if ((previousBinX>0.) ) { mean = (S4_UL + S4_UR + S4_DL + S4_DR)/4.; }
	  else if ((previousBinX==0.) && (nextBinX==0.)) { mean = (previousBinY+S4_UR+S4_DR)/3.; }
	  h_xSec_limit->SetBinContent(ix,iy,mean);

	}

      }
    }
    
  }


  TH2F *h_xSec_limit_5GeV = new TH2F("h_xSec_limit_5GeV","h_xSec_limit_5GeV",405,-12.5,2012.5,405,-12.5,2012.5); h_xSec_limit_5GeV->SetName("h_xSec_limit_5GeV");
  for (UInt_t ix=1; ix<(h_xSec_limit_5GeV->GetNbinsX()+1); ix++) {
    for (UInt_t iy=1; iy<(h_xSec_limit_5GeV->GetNbinsY()+1); iy++) {


      div_t divresult;
      divresultX = div (ix,5);
      divresultY = div (iy,5);
      int binX25GeV = divresultX.quot+1;
      int binY25GeV = divresultY.quot+1;

      h_xSec_limit_5GeV->SetBinContent(ix,iy,h_xSec_limit->GetBinContent(binX25GeV,binY25GeV));

    }
  }


  // legends,etc..
  double xmin = 0.; double xmax = .;
  double ymin = 0.;   double ymax = .;
  TString xaxisName;
  TString yaxisName;

  if (whichSample == "T1tttt") {
    xmin = 400.; xmax = 1400.;
    ymin = 0.;   ymax = 900.;
    xaxisName = "M_{ #tilde{g}} [GeV]";
    yaxisName = "M_{ LSP} [GeV]";
  }

  if (whichSample == "T5tttt") {
    xmin = 875.; xmax = 1400.;
    ymin = 225.; ymax = 1600.;
    xaxisName = "M_{ #tilde{g}} [GeV]";
    yaxisName = "M_{ #tilde{t}} [GeV]";
  }

  if (whichSample == "T1t1t") {
    xmin = 200.; xmax = 800.;
    ymin = 100.;  ymax = 750.;
    xaxisName = "M_{ #tilde{t}} [GeV]";
    yaxisName = "M_{ LSP} [GeV]";
  }

  //  TGraph * getRefXsecGraph(TH2F* limit, string myType, double refMult, int variation=0) 
  TGraph *hLimitCurve      = getRefXsecGraph(h_xSec_limit, whichSample, 1,thUnc);       
  TGraph *hLimitCurve_5GeV = getRefXsecGraph(h_xSec_limit_5GeV, whichSample, 1,thUnc);

  TCanvas* c1 = new TCanvas("c1","c1",700,700);
  c1->SetLogz();
  c1->SetRightMargin(0.196);
  c1->SetLeftMargin(0.135);
  gStyle->SetOptStat(0);

  h_xSec_limit_->GetXaxis()->SetTitle(xaxisName);
  h_xSec_limit_->GetYaxis()->SetTitleOffset(1.4);
  h_xSec_limit_->GetYaxis()->SetTitle(yaxisName);
  h_xSec_limit_->GetZaxis()->SetTitle("95% C.L. upper limit on #sigma (pb) (CL_{s})");
  h_xSec_limit_->GetZaxis()->SetTitleOffset(1.4);
  h_xSec_limit_->GetXaxis()->SetRangeUser(xmin,xmax);
  h_xSec_limit_->GetYaxis()->SetRangeUser(ymin,ymax);
  
  if (whichSample == "T1tttt") { h_xSec_limit_->GetZaxis()->SetRangeUser(0.001,2.2); }
  h_xSec_limit_->Draw("COLZ");

  
  hLimitCurve->SetLineWidth(3);
  hLimitCurve->Draw("same");

  // Xcheck
  TFile *f_SMS1 = TFile::Open("XSEXTION.root","READONLY");
  TH2F *h2DXsec = (TH2F*)f_SMS1->Get("h2DXsec");
  TH2F *h2DXsec_5GeV = (TH2F*)f_SMS1->Get("h2DXsec_5GeV_B");

  TH2F *hExclusion_5GeV = (TH2F*)h_xSec_limit_5GeV->Clone("hExclusion_5GeV");
  TH2F *hExclusion = (TH2F*)h_xSec_limit->Clone("hExclusion");

  double xsecMgluino_1TeV = 0.0243547;
  /*
  for (UInt_t ix=1; ix<(h2DXsec->GetNbinsX()+1); ix++) {
    for (UInt_t iy=1; iy<(h2DXsec->GetNbinsY()+1); iy++) {

      if (h2DXsec->GetBinContent(ix,iy)>0.) { h2DXsec->SetBinContent(ix,iy,xsecMgluino_1TeV); }

    }
  }
  */

  hExclusion->Divide(h_xSec_limit,h2DXsec,1.,1.);
  
  for (UInt_t ix=1; ix<(hExclusion->GetNbinsX()+1); ix++) {
    for (UInt_t iy=1; iy<(hExclusion->GetNbinsY()+1); iy++) {

      if ( (hExclusion->GetBinContent(ix,iy))<=1.) { hExclusion->SetBinContent(ix,iy,1.); }
      else { hExclusion->SetBinContent(ix,iy,0.); }

      if ((h_xSec_limit->GetBinContent(ix,iy))==0.) { hExclusion->SetBinContent(ix,iy,0.); }

    }
  }




  TH2F *hExclusion_5GeV = (TH2F*)h_xSec_limit_5GeV->Clone("hExclusion_5GeV");
  hExclusion_5GeV->Divide(h_xSec_limit_5GeV,h2DXsec_5GeV,1.,1.);

  for (UInt_t ix=1; ix<(hExclusion_5GeV->GetNbinsX()+1); ix++) {
    for (UInt_t iy=1; iy<(hExclusion_5GeV->GetNbinsY()+1); iy++) {

      if ( (hExclusion_5GeV->GetBinContent(ix,iy))<=1.) { hExclusion_5GeV->SetBinContent(ix,iy,1.); }
      else { hExclusion_5GeV->SetBinContent(ix,iy,0.); }

      if ((h_xSec_limit_5GeV->GetBinContent(ix,iy))==0.) { hExclusion_5GeV->SetBinContent(ix,iy,0.); }

    }
  }




  TCanvas *cXTest = new TCanvas("cXTest","cXTest",500,500); 
  gPad->SetGrid();
  hExclusion->GetXaxis()->SetRangeUser(xmin,xmax);
  hExclusion->GetYaxis()->SetRangeUser(ymin,ymax);
  hExclusion->Draw("COLZ");

  TCanvas *cXTest_5GeV = new TCanvas("cXTest_5GeV","cXTest_5GeV",500,500); 
  gPad->SetGrid();
  hExclusion_5GeV->GetXaxis()->SetRangeUser(400.,1400.);
  hExclusion_5GeV->GetYaxis()->SetRangeUser(0.,700.);
  hExclusion_5GeV->Draw("COLZ");

  TCanvas *cXTest_5GeV_TMP = new TCanvas("cXTest_5GeV_TMP","cXTest_5GeV_TMP",700,700); 
  cXTest_5GeV_TMP->SetLogz();
  gPad->SetGrid();
  h_xSec_limit_5GeV->GetXaxis()->SetRangeUser(400.,1500.);
  h_xSec_limit_5GeV->GetYaxis()->SetRangeUser(0.,1125.);
  h_xSec_limit_5GeV->GetZaxis()->SetRangeUser(0.0001,1);
  h_xSec_limit_5GeV->Draw("COLZ");

  hLimitCurve_5GeV->SetLineWidth(3);
  hLimitCurve_5GeV->Draw("same");


  if (saveOutputFile) {

    TString HT, xval;
 
    if (whichHT==500)       { HT = "500";  }
    else if (whichHT==750)  { HT = "750";  }
    else if (whichHT==1000) { HT = "1000"; }

   
    TString rootFileName;
    if (whichSample == "T1tttt") { rootFileName = "SMS_T1tttt_19p4fb_HT500t_MG_Final_ViennaSmoothing.root"; }
    if (whichSample == "T1t1t")  { rootFileName = "SMS_T1t1t_19p4fb_HT500t_MG_Final_Smoothedx5.root"; }
    if (whichSample == "T5tttt") { rootFileName = "SMS_T5tttt_19p4fb_HT500t_MG_Final_Smoothedx5.root"; }

    TFile *fout = new TFile("./rootFiles/"+rootFileName,"UPDATE");
    TString thUnc_ = "";


    if (thUnc==1)       { thUnc_ = "thUp"; }
    else if (thUnc==-1) { thUnc_ = "thDown"; } 


    // save 2D histo
    //    fout->mkdir("Histo2D_XSecUL");
    //    fout->cd("Histo2D_XSecUL");
    h_xSec_limit_->Write("h2D_"+calcWhat+"_"+thUnc_+"_"+HT+"t");
    
    // save limit lines
    //    fout->mkdir("LimitLine");
    //    fout->cd("LimitLine");
    hLimitCurve->Write("LimitLine_"+calcWhat+"_"+thUnc_+"_"+HT+"t");

    // save canvas
    //    fout->mkdir("Canvas");
    //    fout->cd("Canvas");
    c1->Write("Canvas_"+calcWhat+"_"+thUnc_+"_"+HT+"t");

    fout->Close();

    
  }

}
