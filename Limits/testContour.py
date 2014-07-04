#! /usr/bin/env python

import sys, os, commands, math, ROOT, tokenize
from ROOT import TH2F
from array import array


#file = ROOT.TFile("T2tt___sigma_UL_bestexpected.root")
#h    = file.Get("T2tt__expected_xsection_UL") 
file = ROOT.TFile("T2bw_0p25_RL_RR_sigma_UL_bestexpected.root")
h    = file.Get("T2bw_0p25_expected_strength_UL") 
#file = ROOT.TFile("T2tt___sigma_UL_bestexpected.root")
#h    = file.Get("T2tt__expected_strength_UL") 


tg2d = ROOT.TGraph2D(h) # Crea TGraph2D da TH2
contours = [1.0] # Crea array con i valori dei contour che vuoi
a_contours = array('d',contours)
tg2d.GetHistogram().SetContour(1,a_contours) #SetContour(numero contour,array contour)
tg2d.Draw("cont list") # Dummy plotting, serve solo per creare la lista di contour
nb = input('Choose a number: ')
contLevel = tg2d.GetContourList(1.) # Prendi quello che ti interessa
doContour = False
if contLevel.GetSize()>0: doContour = True
#contour = contLevel.First()
print doContour

	
cSmoothed = ROOT.TCanvas("test","test",800,600) # Crea canvas
cSmoothed.cd() 
hrl = cSmoothed.DrawFrame(0.1,0.,1.,0.5); # Commando fondamentale: ti setta il range del Pad, formato (x1,y1,x2,y2) dove (x1,y1) e (x2,y2) sono le coordinate estreme del tuo TH2

h.Draw("colz same") #prima fa il Draw() del TH2 (o TGraph2D)
if doContour:
    for contour in contLevel:
        contour.SetLineWidth(3)
        contour.SetLineStyle(2)
        contour.Draw("L") # poi fa del TGraph del contour, ma senza l'opzione "A"
ROOT.gPad.RedrawAxis();
cSmoothed.Update();                
cSmoothed.SaveAs(".root");
