#! /usr/bin/env python

import ROOT

from ROOT import TH1F

import sys, os, commands, math, copy
from array import array

print "Take the exluded regions and extract the contours"

class Point():
    def __init__(self,x,y):
        self.x = x
        self.y = y

center = Point(540,145)
#center = Point(468,131)

def checkNearbyPoints(histo,i,j):
    if histo.GetBinContent(i+1,j) == 0 or histo.GetBinContent(i-1,j) == 0 or histo.GetBinContent(i,j+1) == 0 or histo.GetBinContent(i,j-1) == 0: 
        return True
    else: return False


def sorterComparator(a,b):
    return ((a.x-center.x)*(b.y-center.y) - (a.y-center.y)*(b.x-center.x)) < 0


def makeComparator(a,b):
    if   sorterComparator(a,b): return 1
    elif sorterComparator(b,a): return -1
    else: return 0


def computeCenter(points):
    sumx = 0
    sumy = 0
    for obj in points:
        sumx = sumx + obj.x
        sumy = sumy + obj.y
    if len(points) == 0:
        return Point(0,0)
    else:
        return Point(sumx/len(points),sumy/len(points))


def addExcludedGraph(model,scenarioX,polschema,variation):
    fBest = ROOT.TFile(model+"_"+scenarioX+"_"+polschema+"_sigma_UL_bestexpected.root","UPDATE")

    typeExlcusion = "_expected_"
    if not variation  == ""    : typeExlcusion = typeExlcusion + variation + "_"
    if variation[0:3] == "obs" : typeExlcusion = "_" + variation  + "_"

    hExcluded = fBest.Get(model+"_"+scenarioX+typeExlcusion+"excludedPoints")
    
    points = []

    for i in range(1,hExcluded.GetNbinsX()+1):
        for j in range(1,hExcluded.GetNbinsY()+1):
            if hExcluded.GetBinContent(i,j) == 1 and checkNearbyPoints(hExcluded,i,j):
                points.append(Point(hExcluded.GetXaxis().GetBinCenter(i),hExcluded.GetYaxis().GetBinCenter(j)))

    center = computeCenter(points)
    print model, scenarioX, polschema, variation, "center: {0:.1f} {1:.1f}".format(center.x, center.y)

    sortedPoints = sorted(points, cmp=makeComparator)
    mstops = []
    mlsps  = []
    for obj in sortedPoints:
        mstops.append(obj.x)
        mlsps.append(obj.y)
    graph = ROOT.TGraph(len(mstops),array('d',mstops),array('d',mlsps))
    #hExcluded.Draw("colz")
    #graph.Draw("PL*")
    #input()
    fBest.cd()
    graph.Write(model+"_"+scenarioX+typeExlcusion+"excludedRegion")
    fBest.Close()




from modelclass import *
models = [Model("T2tt", [""]                  , [""]),
          Model("T2bw", ["0p25","0p50","0p75"], ["LL_LR"])]#,"RL_RR"])]

variations = ["","p1s","m1s","observed","observed_p1s","observed_m1s"]
#variations = ["","p1s","m1s","L","R","observed","observed_p1s","observed_m1s"]

for model in models:
    for scenarioX in model.scenariosX:
        for polschema in model.regions:
            for variation in variations:
                addExcludedGraph(model.model,scenarioX,polschema,variation)
    
