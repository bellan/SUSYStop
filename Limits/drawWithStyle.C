#include <TROOT.h>
#include <TMultiGraph.h>
#include "setTDRStyle.C"

void paletteColdToHot(TStyle *style, TString type = "TChiWX");
void embellish(TH2F*, const TString& titleXaxis, const TString& titleYaxis, const TString& titleZaxis="");
void drawCMSInscription();

void paletteColdToHot(TStyle *style, TString type){
  // ********************** //
  // Set Display properties
  // ********************** //
  
  //TStyle *tdrStyle = setTDRStyle();
  //tdrStyle->SetCanvasDefW(700);
  //paletteColdToHot(tdrStyle);

  // Sue Ann' style
  if(type == "SAK"){
    const Int_t     NRGBs         = 6;
    const Int_t     NCont         = 255;
    Double_t        stops [NRGBs] = { 0.00, 0.25, 0.5, 0.7  , 0.85, 1.00 };
    Double_t        red   [NRGBs] = { 0.00, 0.05, 0.65, 1.00, 1.00, 1.00 };
    Double_t        green [NRGBs] = { 0.00, 0.85, 1.00, 1.00, 0.75, 0.00 };
    Double_t        blue  [NRGBs] = { 0.51, 1.00, 0.15, 0.30, 0.20, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    style->SetNumberContours(NCont);
  }
  
  // Fedor' style
  if(type == "TChiWX"){
    const Int_t     NRGBs         = 5;
    const Int_t     NCont         = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.50, 0.50, 1.00, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.50, 1.00, 1.00, 0.60, 0.50 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.50, 0.40, 0.50 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    style->SetNumberContours(NCont);
  }

  if(type == "TChiWX2"){
    const Int_t     NRGBs         = 5;
    const Int_t     NCont         = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 1.00, 0.50, 0.40, 0.50 };
    Double_t green[NRGBs] = { 0.50, 1.00, 1.00, 0.60, 0.50 };
    Double_t blue[NRGBs]  = { 0.50, 0.50, 1.00, 1.00, 1.00 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    style->SetNumberContours(NCont);
  }



  if(type == "TChiWX3"){
    const Int_t     NRGBs         = 9;
    const Int_t     NCont         = 255;
    Double_t stops[NRGBs] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.75, 0.750, 0.625, 0.500, 0.5, 0.50, 0.5, 0.5    };
    Double_t green[NRGBs] = { 0.25, 0.250, 0.750, 1.000, 1.000, 0.75, 0.50, 0.50, 0.50 };
    Double_t blue[NRGBs]  = { 0.  , 0.1, 0.250, 0.365, 0.500, 0.625, 0.75, 0.8, 1.00 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    style->SetNumberContours(NCont);
  }



  if(type == "Test"){
    const Int_t     NRGBs         = 5;
    const Int_t     NCont         = 255;
    //Double_t stops[NRGBs] = { 1.00, 0.84, 0.61, 0.34, 0.00 };
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 1.00, 1.00, 0.50, 0.50 };
    Double_t green[NRGBs] = { 0.50, 0.60, 1.00, 1.00, 0.50 };
    Double_t blue[NRGBs]  = { 0.50, 0.40, 0.50, 1.00, 1.00 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    style->SetNumberContours(NCont);
  }

  //gStyle->SetPaintTextFormat(".2f");
}


  //paletteGrayGreen(10,20);

//_____________________________________________________________________________
void paletteGrayGreen(const Int_t numGray, const Int_t numGreen)
{
  static const Float_t    GRAY_HIGH   = 0.9f;
  static const Float_t    GRAY_LOW    = 0.1f;
  static const Float_t    RED_HIGH    = 1.f;
  static const Float_t    RED_LOW     = 0.44f;
  static const Float_t    GREEN_HIGH  = 0.9f;
  static const Float_t    GREEN_LOW   = 0.5f;
  static const Float_t    BLUE_HIGH   = 0.3f;
  static const Float_t    BLUE_LOW    = 0.32f;

  const Float_t           grayStep    = (GRAY_HIGH  - GRAY_LOW ) / (numGray  - 1);
  const Float_t           redStep     = (RED_HIGH   - RED_LOW  ) / (numGreen - 1);
  const Float_t           greenStep   = (GREEN_HIGH - GREEN_LOW) / (numGreen - 1);
  const Float_t           blueStep    = (BLUE_HIGH  - BLUE_LOW ) / (numGreen - 1);
  const Int_t             numLevels   = numGreen + numGray;
  std::vector<Int_t>      colors(numLevels);

  // May need to initialize the color system
  TSeqCollection*         allColors   = gROOT->GetListOfColors();
  if (!allColors || allColors->GetEntries() == 0) {
    TColor::InitializeColors();
    allColors   = gROOT->GetListOfColors();
    if (allColors == 0 || allColors->GetEntries() == 0)
      cout << "Presentation::paletteInverseGrayscale():  Something is wrong with the color system!" << endl;
  }


  Int_t                   colorIndex  = TMath::Max(1000, 1 + dynamic_cast<TColor*>(allColors->Last())->GetNumber());
  Int_t                   iColor      = 0;
  
  for (Float_t value = GRAY_HIGH; iColor < numGray; ++iColor, ++colorIndex, value -= grayStep) {
    new TColor(colorIndex, value, value, value);
    colors[iColor]        = colorIndex;
  } // end loop over generated colors

  for (Float_t red = RED_LOW, blue = BLUE_LOW, green = GREEN_LOW; iColor < numLevels; ++iColor, ++colorIndex, red += redStep, green += greenStep, blue += blueStep) {
    new TColor(colorIndex, red, green, blue);
    colors[iColor]        = colorIndex;
  } // end loop over generated colors


  TColor::SetPalette(numLevels, &(colors[0]));
  gStyle->SetNumberContours(numLevels);
}




void embellish(TH2F* histo, const TString& titleXaxis, const TString& titleYaxis, const TString& titleZaxis){
  
 
  histo->GetXaxis()->SetTitle(titleXaxis);
  histo->GetYaxis()->SetTitle(titleYaxis);
  histo->GetZaxis()->SetTitle(titleZaxis);

  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetTitleOffset(1.5);
  histo->GetZaxis()->SetTitleOffset(1.5);
  
  histo->SetStats(kFALSE);
}

void setStyle(TMultiGraph* outline, int color, int style, int width = 5){
  TList *l = outline->GetListOfGraphs();
  TIter next(l);
  while (TObject *obj = next()){
    ((TGraph*)obj)->SetLineWidth(width);
    ((TGraph*)obj)->SetLineColor(color);
    ((TGraph*)obj)->SetLineStyle(style);
  }
}

void setStyle(TGraph2D* obj, int color, int style, int width = 5){
  obj->SetLineWidth(width);
  obj->SetLineColor(color);
  obj->SetLineStyle(style);
}

void setStyle(TGraph* obj, int color, int style, int width = 5){
  obj->SetLineWidth(width);
  obj->SetLineColor(color);
  obj->SetLineStyle(style);
}



void drawCMSInscription(){
  
  TLatex l;
  l.SetTextAlign(12);
  l.SetTextSize(0.035);
  l.SetTextFont(132);
  
  l.SetNDC();

  l.DrawLatex(0.20,0.972,"CMS Preliminary");
  l.DrawLatex(0.17,0.935,"#sqrt{s} = 8 TeV   L = 19.03 fb^{-1}");
}
