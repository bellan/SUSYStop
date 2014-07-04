#! /usr/bin/env python

import sys, os, commands, math, ROOT, tokenize, collections
from ROOT import TH2F



def getValues(model, region, limitType):

    inputdir  = "CLs_"+limitType
    model_region = model+"_"+region
    failure1, output1 = commands.getstatusoutput("ls {0:s}/*{1:s}*.root".format(inputdir,model_region))
    rootfiles = output1.split()
    values = {}
    print model_region, len(rootfiles)
    for file in rootfiles:
        file = ROOT.TFile(file)
        tree = file.Get("limit")
        entries = tree.GetEntries()
        for jentry in xrange(entries):           
            # get the next tree in the chain and verify
            ientry = tree.LoadTree(jentry)
            if ientry < 0: break
            # copy next entry into memory and verify
            nb = tree.GetEntry(jentry)
            if nb<=0: continue
            # use the values directly from the tree
            #    print mychain.MHT, weight
            values[round(tree.mh,2)] = round(tree.limit,4) 
    values = collections.OrderedDict(sorted(values.items()))
    return values





def fillPlots(model,scenario,region,listunc):
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

    
    observed     = getValues(model, region, "observed")
    expected     = getValues(model, region, "expected_0.50")
    expected_p1s = getValues(model, region, "expected_0.84")
    expected_m1s = getValues(model, region, "expected_0.16")

    for i in range(1,hTheounc.GetNbinsX()+1):
        for j in range(1,hTheounc.GetNbinsY()+1):
            if hTheounc.GetBinContent(i,j) == 0: continue
            unc = hTheounc.GetBinContent(i,j)
            for l in listunc:
                if l > unc:
                    hObs.SetBinContent(i,j,observed[l+1])
                    hExp.SetBinContent(i,j,expected[l+1])
                    hExp1p.SetBinContent(i,j,expected_p1s[l+1])
                    hExp1m.SetBinContent(i,j,expected_m1s[l+1])
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



#for model in models:
#    for region in model.regions:
#        values = getValues(model,region,limitType)
#        print values[1.15]


from modelclass import *

models = [Model("T2bw", ["0p25","0p50","0p75"], ["LX","LM","MXHM","VHM","HXHM"]),
          Model("T2tt", [""]                  , ["LM","MM","HM","VHM"])]

for model in models:
    for scenarioX in model.scenariosX:
        for region in model.regions:
            #print model.model, scenarioX, region
            fillPlots(model.model, scenarioX, region, [x/100. for x in range(0,125,5)])

