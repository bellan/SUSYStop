
{
  TString model = "T2bw";
  TString scenarioX = "0p25";

  TFile *fn = new TFile("bestexpected_"+model+"_"+scenarioX+"_stain.root");
  TFile *fl = new TFile("bestexpected_L_"+model+"_"+scenarioX+"_stain.root");
  TFile *fr = new TFile("bestexpected_R_"+model+"_"+scenarioX+"_stain.root");

  TH2F* hn = fn->Get("excludedPoints");
  TH2F* hl = fl->Get("excludedPoints");
  TH2F* hr = fr->Get("excludedPoints");

  hn->SetFillColor(4);
  hn->Draw("box");

  hl->SetFillColor(2);
  hl->Draw("boxsame");

  hr->SetFillColor(3);
  hr->Draw("boxsame");
}
