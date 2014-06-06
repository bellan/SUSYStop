#! /usr/bin/env python

import ROOT

from ROOT import TH1F

import sys, os, commands, math, copy

model = sys.argv[1]
scenarioX = ""
pol = ""

if model == "T2bw":
    scenarioX = sys.argv[2]
    pol = "LL_LR"

fn = ROOT.TFile("bestexpected_"+model+"_"+scenarioX+"_"+pol+"_stain.root")
fl = ROOT.TFile("bestexpected_L_"+model+"_"+scenarioX+"_"+pol+"_stain.root")
fr = ROOT.TFile("bestexpected_R_"+model+"_"+scenarioX+"_"+pol+"_stain.root")

hn = fn.Get("excludedPoints")
hl = fl.Get("excludedPoints")
hr = fr.Get("excludedPoints")

canvas = ROOT.TCanvas("Stain_"+model,"Stain_"+model)
canvas.cd()

#hn.SetFillColor(4)
#hn.SetFillStyle(3145)
#hn.Draw("box")

#hl.SetFillStyle(3154)
#hl.Draw("boxsame")
hl.SetContour(2)
hl.SetLineColor(ROOT.kSpring)
hl.SetLineWidth(5)
hl.Draw("cont3")


#hr.SetFillColor(3)
#hr.SetFillStyle(3150)
#hr.Draw("boxsame")
hr.SetContour(2)
hr.SetLineColor(ROOT.kOrange)
hr.SetLineWidth(5)
hr.Draw("cont3same")

if model == "T2bw":
    pol = "RL_RR"
    frl = ROOT.TFile("bestexpected_L_"+model+"_"+scenarioX+"_"+pol+"_stain.root")
    frr = ROOT.TFile("bestexpected_R_"+model+"_"+scenarioX+"_"+pol+"_stain.root")
    hrl = frl.Get("excludedPoints")
    hrr = frr.Get("excludedPoints")
    hrl.SetContour(2)
    hrl.SetLineColor(ROOT.kBlack)
    hrl.SetLineWidth(5)
    hrl.Draw("cont3same")

    hrr.SetContour(2)
    hrr.SetLineColor(ROOT.kBlue)
    hrr.SetLineWidth(5)
    hrr.Draw("cont3same")

hn.SetContour(2)
hn.SetLineColor(ROOT.kRed)
hn.SetLineWidth(5)
hn.Draw("cont3same")


#input()

out = ROOT.TFile(model+"_"+scenarioX+"_metro.root","RECREATE")
out.cd()
canvas.Write(model+"_"+scenarioX)
out.Close()




