#! /usr/bin/env python

import sys, os, commands, math, ROOT, tokenize
from ROOT import TH2F


limitType = "observed"
inputdir  = "CLs_"+limitType

from modelclass import *

models = [Model("T2bw", [""], ["LX","LM","MXHM","VHM","HXHM"]),
          Model("T2tt", [""], ["LM","MM","HM","VHM"])]

for model in models:
    for region in model.regions:
        model_region = model.model+"_"+region
        print model_region
        failure1, output1 = commands.getstatusoutput("ls {0:s}/*{1:s}*.root".format(inputdir,model_region))
        rootfiles = output1.split()
        command = "hadd {0:s}/{0:s}_{1:s}.root ".format(inputdir, model_region)
        for file in rootfiles:
            command = command + " " + file
        print command
        failure2, output2 = commands.getstatusoutput(command)
