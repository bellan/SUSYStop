#! /usr/bin/env python

import ROOT

from ROOT import TH1F

import sys, os, commands, math, copy

print "Produce sigma UL for direct stop search"

def calculatedSigma_UL(model, scenarioX, region, polschema):

    #print "--------------------------------------"
    print model, scenarioX, region, polschema
    #print "--------------------------------------"

    lumi = 19030

    hEff = ROOT.TH1F("FakeEff","FakeEff",1,0,1)
    
    fname = ""
    if model == "T2tt": 
        fname = "theouncertainty/"+model+"_"+region+"_cvAndSys.root"
    if model == "T2bw": 
        fname = "theouncertainty/"+model+"_"+scenarioX+"_"+region+"_cvAndSys.root"

    f = ROOT.TFile(fname)
    hEff  = f.Get("Efficiency")

    hEff.GetZaxis().SetRangeUser(0,1)

    # polarization inputs
    #fPolName = "Polarization/pol_"+model+"_"+scenarioX+"_"+polschema+"_Had.root"
    fPolName = "polarization/pol_"+model+"_"+scenarioX+"_"+polschema+".root"

    fPol = ROOT.TFile(fPolName)
    hL_ratio = fPol.Get(model.upper()+"_"+region+"__VARIATIONAXIS_Down_ratio") # FIXME 
    hR_ratio = fPol.Get(model.upper()+"_"+region+"__VARIATIONAXIS_Up_ratio") # FIXME

    hExp = ROOT.TH2F(model+"_expected_xsection_UL",
                     model+"_expected_xsection_UL",
                     hEff.GetXaxis().GetNbins(), hEff.GetXaxis().GetXmin(), hEff.GetXaxis().GetXmax(),
                     hEff.GetYaxis().GetNbins(), hEff.GetYaxis().GetXmin(), hEff.GetYaxis().GetXmax()) 
 
    hExp_p1s = hExp.Clone(model+"_expected_p1s_xsection_UL")
    hExp_m1s = hExp.Clone(model+"_expected_m1s_xsection_UL")
    hObs     = hExp.Clone(model+"_observed_xsection_UL")

    # polarizations output
    hExpL     = hExp.Clone(model+"_expected_L_xsection_UL")
    hExpR     = hExp.Clone(model+"_expected_R_xsection_UL")
    hObsL     = hExp.Clone(model+"_observed_L_xsection_UL")
    hObsR     = hExp.Clone(model+"_observed_R_xsection_UL")

    # input plots representing the UL on the events
    fEvULname = model+"_"+scenarioX+"_"+region+"_nevents_UL.root"
    fEvUL    = ROOT.TFile(fEvULname)

    hObsEv   = fEvUL.Get("ULObservedEvents")
    hExpEv   = fEvUL.Get("ULExpectedEvents")
    hExp1pEv = fEvUL.Get("ULExpected1pEvents")
    hExp1mEv = fEvUL.Get("ULExpected1mEvents")

    fout = ROOT.TFile(str(model)+"_"+str(scenarioX)+"_"+str(region)+"_sigma_UL.root","RECREATE")

    for i in range(1,hEff.GetNbinsX()+1):
        for j in range(1,hEff.GetNbinsY()+1):
            eff = hEff.GetBinContent(i,j)

            if not eff == 0:
                exp_NUL      =  hExpEv.GetBinContent(i,j)
                exp_NUL_p1s  =  hExp1pEv.GetBinContent(i,j)
                exp_NUL_m1s  =  hExp1mEv.GetBinContent(i,j)
                obs_NUL      =  hObsEv.GetBinContent(i,j)

                #print "Bin [{0:.0f},{1:.0f}]: exp UL = {2:.3f}, eff = {3:.3f}, exp sigma UL = {4:.3f}".format(hEff.GetXaxis().GetBinCenter(i), hEff.GetYaxis().GetBinCenter(j), exp_NUL, eff, exp_NUL/(eff*lumi))
                hExp.SetBinContent(i,j,exp_NUL/(eff*lumi))
                hExp_p1s.SetBinContent(i,j,exp_NUL_p1s/(eff*lumi))
                hExp_m1s.SetBinContent(i,j,exp_NUL_m1s/(eff*lumi))
                hObs.SetBinContent(i,j,obs_NUL/(eff*lumi))

                # FIXME!!! Temporary bug fix for  LSP mass == 0 GeV
                #if model == "T2bw" and hExp.GetYaxis().GetBinCenter(j) == 0:
                #    hExp.SetBinContent(i,j,0)
                #    hExp_p1s.SetBinContent(i,j,0)
                #    hExp_m1s.SetBinContent(i,j,0)
                #    hObs.SetBinContent(i,j,0)


                #polarizarions
                # polarization files have a different bin range. x-axis starts from 0 instead of 200 GeV as for the other histos
                i_pol = hL_ratio.GetXaxis().FindBin(hExp.GetXaxis().GetBinCenter(i))
                j_pol = hL_ratio.GetYaxis().FindBin(hExp.GetYaxis().GetBinCenter(j))

                if not hL_ratio.GetBinContent(i_pol,j_pol) == 0:
                    hExpL.SetBinContent(i,j,hExp.GetBinContent(i,j)/hL_ratio.GetBinContent(i_pol,j_pol))
                    hObsL.SetBinContent(i,j,hObs.GetBinContent(i,j)/hL_ratio.GetBinContent(i_pol,j_pol))
                else:
                    hExpL.SetBinContent(i,j,0)
                    hObsL.SetBinContent(i,j,0)

                if not hR_ratio.GetBinContent(i_pol,j_pol) == 0:
                    hExpR.SetBinContent(i,j,hExp.GetBinContent(i,j)/hR_ratio.GetBinContent(i_pol,j_pol))
                    hObsR.SetBinContent(i,j,hObs.GetBinContent(i,j)/hR_ratio.GetBinContent(i_pol,j_pol))
                else:
                    hExpR.SetBinContent(i,j,0)
                    hObsR.SetBinContent(i,j,0)
               
            else:
                hExp.SetBinContent(i,j,0)
                hExp_p1s.SetBinContent(i,j,0)
                hExp_m1s.SetBinContent(i,j,0)
                hObs.SetBinContent(i,j,0)

                hExpL.SetBinContent(i,j,0)
                hExpR.SetBinContent(i,j,0)
                hObsL.SetBinContent(i,j,0)
                hObsR.SetBinContent(i,j,0)

    fout.cd()
    hEff.Write()
    hExp.Write()
    hExp_m1s.Write()
    hExp_p1s.Write()
    hObs.Write()

    hExpL.Write()
    hExpR.Write()
    hObsL.Write()
    hObsR.Write()

    fout.Close()


model = sys.argv[1]
models = []
models.append(model)

scenariosX = []
regions = []
polschema = ""

if model == "T2bw":
    scenariosX = ["0p25", "0p50", "0p75"]
    regions    = ["LX","LM","MXHM","VHM","HXHM"]
    #pol = "LL_LR" or "LR_RR"
    polschema = sys.argv[2]

if model == "T2tt":
    scenariosX = [""]
    regions    = ["LM","MM","HM","VHM"]

#print model, "polarization schema: ", polschema

for model in models:
    for scenarioX in scenariosX:
        for region in regions:
            calculatedSigma_UL(model, scenarioX, region, polschema)


# Take the best out of the three regions
print "Look for the best, point by point, for all regions. Scan each region."
for model in models:
    for scenarioX in scenariosX:
        files  = []
        histos = []
        histos1p = []
        histos1m = []

        histosL = []
        histosR = []
        
        for region in regions:
            f = ROOT.TFile(str(model)+"_"+str(scenarioX)+"_"+str(region)+"_sigma_UL.root")
            #print str(model),str(scenarioX),str(region)
            
            histos.append(copy.deepcopy(  f.Get(model+"_expected_xsection_UL")))
            histos1p.append(copy.deepcopy(f.Get(model+"_expected_p1s_xsection_UL")))
            histos1m.append(copy.deepcopy(f.Get(model+"_expected_m1s_xsection_UL")))

            histosL.append(copy.deepcopy(f.Get(model+"_expected_L_xsection_UL")))
            histosR.append(copy.deepcopy(f.Get(model+"_expected_R_xsection_UL")))

        best = ROOT.TH2F(str(model)+"_"+str(scenarioX)+"_expected_xsection_UL",
                         str(model)+"_"+str(scenarioX)+"_expected_xsection_UL",
                         histos[0].GetXaxis().GetNbins(), histos[0].GetXaxis().GetXmin(), histos[0].GetXaxis().GetXmax(),
                         histos[0].GetYaxis().GetNbins(), histos[0].GetYaxis().GetXmin(), histos[0].GetYaxis().GetXmax())
        
        best1p = best.Clone(str(model)+"_"+str(scenarioX)+"_expected_p1s_xsection_UL")
        best1p.SetTitle    (str(model)+"_"+str(scenarioX)+"_expected_p1s_xsection_UL")
        best1m = best.Clone(str(model)+"_"+str(scenarioX)+"_expected_m1s_xsection_UL")
        best1m.SetTitle    (str(model)+"_"+str(scenarioX)+"_expected_m1s_xsection_UL")

        obs = best.Clone(str(model)+"_"+str(scenarioX)+"_observed_xsection_UL")
        obs.SetTitle    (str(model)+"_"+str(scenarioX)+"_observed_xsection_UL")

        #polarizations
        bestL = best.Clone(str(model)+"_"+str(scenarioX)+"_expected_L_xsection_UL")
        bestL.SetTitle    (str(model)+"_"+str(scenarioX)+"_expected_L_xsection_UL")
        bestR = best.Clone(str(model)+"_"+str(scenarioX)+"_expected_R_xsection_UL")    
        bestR.SetTitle    (str(model)+"_"+str(scenarioX)+"_expected_R_xsection_UL")    

        bestregionP = best.Clone(str(model)+"_"+str(scenarioX)+"_bestRegion")
        bestregionP.SetTitle    (str(model)+"_"+str(scenarioX)+"_bestRegion")

        for x in range(1,best.GetNbinsX()+1):
            for y in range(1,best.GetNbinsY()+1):
                bestpoint = 10000
                bestregion = -1
                for index,histo in enumerate(histos):
                    if histo.GetBinContent(x,y) < bestpoint and not histo.GetBinContent(x,y) == 0:
                        bestpoint = histo.GetBinContent(x,y)
                        bestregion = index

                if bestregion == -1 and not histo.GetBinContent(x,y) == 0: print "WARNING region not recognized!!"
                best.SetBinContent(x,y,histos[bestregion].GetBinContent(x,y))
                best1p.SetBinContent(x,y,histos1p[bestregion].GetBinContent(x,y))
                best1m.SetBinContent(x,y,histos1m[bestregion].GetBinContent(x,y))
        
                bestL.SetBinContent(x,y,histosL[bestregion].GetBinContent(x,y))
                bestR.SetBinContent(x,y,histosR[bestregion].GetBinContent(x,y))

                if not bestpoint == 0: bestregionP.SetBinContent(x,y,bestregion+1)
        fout = ROOT.TFile(str(model)+"_"+str(scenarioX)+"_"+polschema+"_sigma_UL_bestexpected.root","RECREATE")

                    

        fout.cd()
        best.Write()
        best1p.Write()
        best1m.Write()
        obs.Write()
        bestregionP.Write()

        bestL.Write()
        bestR.Write()

        fout.Close()
