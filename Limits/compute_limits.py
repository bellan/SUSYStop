#! /usr/bin/env python

import ROOT

from ROOT import TH1F

import sys, os, commands, math, copy, ast


### ---------------------------------------------------------------------------------------------------------------- ###
print "Prepare inputs for limits calculation"
### ---------------------------------------------------------------------------------------------------------------- ###

os.popen('rm T2*sigma_UL_bestexpected.root')
intResults = "intermediateResults"
if os.path.exists(intResults):
    os.popen('rm -r "%s"' %intResults)
os.popen('mkdir "%s"' %intResults)

failure = ""
limitType = sys.argv[1] if len(sys.argv) >1 else "expected"
debug = ast.literal_eval(sys.argv[2]) if len(sys.argv) == 3 else False


### ---------------------------------------------------------------------------------------------------------------- ###

failure1, output1 = commands.getstatusoutput("./compute_nevents_UL_fromROOTfiles.py")
if debug: print output1
failure = str(failure1)

### ---------------------------------------------------------------------------------------------------------------- ###

from modelclass import *
models = [Model("T2tt", [""]                  , [""]),
          Model("T2bw", ["0p25","0p50","0p75"], ["LL_LR"])]#,"RL_RR"])]         

for model in models:
    for polschema in model.regions:
        failure2, output2 = commands.getstatusoutput("./compute_sigma_UL.py {0:s} {1:s}".format(model.model,polschema))
        if debug: print output2
        failure = failure + " " + str(failure2)

print "Move intermediate results in",intResults
os.popen('mv T2*UL.root "%s"' %intResults)

### ---------------------------------------------------------------------------------------------------------------- ###
print "Produces limits"
### ---------------------------------------------------------------------------------------------------------------- ###

root = "root -l -e -q "
for model in models:
    for scenarioX in model.scenariosX:
        for pol in model.regions: # polarizations in this case!
            macro = "'addExcludedRegions.C+("+'"{0:s}","{1:s}","{2:s}")'.format(model.model,scenarioX,pol)+"'"
            command = root+macro
            print command
            failure3, output3 = commands.getstatusoutput(command)
            failure = failure + " " + str(failure3)
            if debug: print output3

### ---------------------------------------------------------------------------------------------------------------- ###

failure4, output4 = commands.getstatusoutput("./addExcludedGraphs.py")
print "./addExcludedGraphs.py"
if debug: print output4
failure = failure + " " + str(failure4)

### ---------------------------------------------------------------------------------------------------------------- ###

#doPolChoices = [True,False]
doPolChoices = [False]
preliminary = True

for model in models:
    for scenarioX in model.scenariosX:
        for polschema in model.regions: # polarizations in this case!
            for doPol in doPolChoices: 
                macro = "'makePlots_smoothing.C+("+'"{0:s}","{1:s}","{2:s}",{3:b},"{4:s}",{5:b})'.format(model.model,scenarioX,polschema,doPol,limitType,preliminary)+"'"
                command = root+macro
                print command
                failure4, output4 = commands.getstatusoutput(command)
                failure = failure + " " + str(failure4)
                #if debug: print output4
                print output4

### ---------------------------------------------------------------------------------------------------------------- ###

print "Exit status codes:", failure
