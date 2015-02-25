void testContour()
{
  TFile *file = new TFile("T2tt___sigma_UL_bestexpected.root");
  TH2F  *h    = (TH2F*)file->Get("T2tt__expected_strength_UL");
  
  TGraph2D *tg2d = new TGraph2D(h); //Crea TGraph2D da TH2
  Double_t contours[1]; // Crea array con i valori dei contour che vuoi
  contours[0] = 1.0;
  tg2d->GetHistogram()->SetContour(1,contours);  //SetContour(numero contour,array contour)
  //tg2d->Draw("cont list"); //Dummy plotting, serve solo per creare la lista di contour
  TList *contLevel = tg2d->GetContourList(1.); // Prendi quello che ti interessa
  
  if(contLevel->GetSize()>0) doContour = true;
  cout<< doContour << endl;
  
  TCanvas *cSmoothed = new TCanvas("test","test",800,600); // Crea canvas
  cSmoothed->cd();
  //TH1F* hrl = cSmoothed->DrawFrame(0.1,0.,1.,0.5); // Commando fondamentale: ti setta il range del Pad, formato (x1,y1,x2,y2) dove (x1,y1) e (x2,y2) sono le coordinate estreme del tuo TH2
  
  h->Draw("colz same"); //prima fa il Draw() del TH2 (o TGraph2D)
  if (doContour){
    
    TIter next(contLevel);
    TObject *contour = 0;
    while (contour = next()){
      ((TGraph2D*)contour)->SetLineWidth(3);
      ((TGraph2D*)contour)->SetLineStyle(2);
      ((TGraph2D*)contour)->Draw("L"); // poi fa del TGraph del contour, ma senza l'opzione "A"
    }
  }
  //ROOT::gPad::RedrawAxis();
  cSmoothed->Update();                
  cSmoothed->SaveAs(".root");
}
