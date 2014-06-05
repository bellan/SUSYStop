#! /usr/bin/env python


import sys, os, commands, math, ROOT, tokenize
from ROOT import TH2F


def getValues(filename, name, printTable = False):
    file = open(filename,'r')
    array = {'Obs': -1,'Exp': -1,'Expm1s': -1,'Expp1s': -1}
    for line in file:
        # find the line with the input

        found =  line.find("Observed")
        if found >= 0:
            # split the line using comma as separator
            spline =  line.split(" ")
            obs =  float(spline[len(spline)-1])
            # print "observed", obs
            array["Obs"] = obs

        found =  line.find("Expected 50.0")
        if found >= 0:
            spline =  line.split(" ")
            exp    =  float(spline[len(spline)-1])
            # print "expected", exp
            array["Exp"] = exp

        found =  line.find("Expected 16.0")
        if found >= 0:
            spline =  line.split(" ")
            exp_m1s    =  float(spline[len(spline)-1])
            # print "-1s expected",exp_m1s
            array["Expm1s"] = exp_m1s

        found =  line.find("Expected 84.0")
        if found >= 0:
            spline =  line.split(" ")
            exp_p1s    =  float(spline[len(spline)-1])
            # print "+1s expected",exp_m1s
            array["Expp1s"] = exp_p1s
        
        found = line.find("@ 95% CL")
        if found >= 0:
            spline  =  line.split(" ")
            exp     =  float(spline[3])
            exp_p1s =  exp + float(spline[5])
            exp_m1s =  exp - float(spline[5])
            array["Exp"] = exp
            array["Expp1s"] = exp_p1s
            array["Expm1s"] = exp_m1s
            
    if printTable: 
        print "{0:s}:     & {1:.1f} & ${2:.1f}_".format(name,array["Obs"],array["Exp"])+"{"+"{0:.1f}".format(array["Expm1s"]-array["Exp"])+"}^{"+"+ {0:.1f}".format(array["Expp1s"]-array["Exp"])+"}$ \\"+"{0:s}".format("\\")
    return array





def fillPlots(model,scenario,region,listunc,printTable=False):
    fTheouncName = ""
    if model == "T2bw":
        fTheouncName = "theouncertainty/"+model+"_"+scenario+"_"+region+"_cvAndSys.root"
    if model == "T2tt":
        fTheouncName = "theouncertainty/"+model+"_"+region+"_cvAndSys.root"

    fTheounc = ROOT.TFile(fTheouncName)
    hTheounc = fTheounc.Get("TotalSysUnc")
    
    hTheouncCoarse = ROOT.TH2F("TheoUncertainty","Theoretical uncertainty",hTheounc.GetXaxis().GetNbins(),hTheounc.GetXaxis().GetXmin(), hTheounc.GetXaxis().GetXmax(),hTheounc.GetYaxis().GetNbins(), hTheounc.GetYaxis().GetXmin(), hTheounc.GetYaxis().GetXmax())


    hObs = hTheouncCoarse.Clone("ULObservedEvents")
    hObs.SetTitle("UL observed events")

    hExp = hTheouncCoarse.Clone("ULExpectedEvents")
    hExp.SetTitle("UL expected events")

    hExp1p = hTheouncCoarse.Clone("ULExpected1pEvents")
    hExp1p.SetTitle("UL expected events (+1s)")

    hExp1m = hTheouncCoarse.Clone("ULExpected1mEvents")
    hExp1m.SetTitle("UL expected events (-1s)")

    for i in range(1,hTheounc.GetNbinsX()+1):
        for j in range(1,hTheounc.GetNbinsY()+1):
            if hTheounc.GetBinContent(i,j) == 0: continue
            unc = hTheounc.GetBinContent(i,j)
            for l in listunc:
                if l > unc:
                    rname = ""
                    values = getValues("results/"+model+"_"+region+"_"+str(1+l)+".txt",region)
                    if printTable: print values
                    hObs.SetBinContent(i,j,values["Obs"])
                    hExp.SetBinContent(i,j,values["Exp"])
                    hExp1p.SetBinContent(i,j,values["Expp1s"])
                    hExp1m.SetBinContent(i,j,values["Expm1s"])
                    hTheouncCoarse.SetBinContent(i,j,l)
                    break

    fOutname = model+"_"+scenarioX+"_"+region+"_nevents_UL.root"
    
    fOut = ROOT.TFile(fOutname,"RECREATE")
    fOut.cd()
    hObs.Write()
    hExp.Write()
    hExp1m.Write()
    hExp1p.Write()
    hTheouncCoarse.Write()
    fOut.Close()



from modelclass import *

models = [Model("T2bw", ["0p25","0p50","0p75"], ["LX","LM","MXHM","VHM","HXHM"]),
          Model("T2tt", [""]                  , ["LM","MM","HM","VHM"])]

for model in models:
    for scenarioX in model.scenariosX:
        for region in model.regions:
            #print model.model, scenarioX, region
            fillPlots(model.model, scenarioX, region, [x/100. for x in range(0,125,5)])


