#! /usr/bin/env python

import ROOT

from ROOT import TH1F

import sys, os, commands, math, ast

def scan(file, blind):

    failure, output = commands.getstatusoutput('basename "%s"' %file)
    template = output

    # base name for the datacard based on the template
    inputbase =  "scenarios/"+template[0:len(template)-4]
    # basename for the output of the limit calculation
    outputbase = "results/"+template[17:len(template)-4]

    # Theo unc scenarios
    list = [1+x/100. for x in range(0,125,5)]
    
    print "Scanning ", template, "using", list

    
    for i in list:
        
        # prepare the specific data card
        input = inputbase+"_"+str(i)+".txt"
        sed = "sed 's,1.XX,{0:.2f},g'  {1:s} > {2:s}".format(i,file,input)
        commands.getstatusoutput(sed)
        
        output = outputbase+"_"+str(i)+".txt"
        blindOption = "--run blind" if blind else ""
        #run_limits =  "combine -M HybridNew --rule CLs --testStat=LHC --generateNuis=0 --fitNuis=0 --generateExt=1 {0:s} > {1:s}".format(input,output)
        run_limits = "combine -M Asymptotic --rMax 100 {0:s} {1:s} > {2:s}".format(blindOption, input, output)
        print run_limits
        commands.getstatusoutput(run_limits)


# Run blind if nothing is specified, otherwise produce observed limits too
blind = ast.literal_eval(sys.argv[1]) if len(sys.argv) == 2 else True


# Get the templates
ls = "ls datacards/limits_datacards_*.txt"
failure,output =  commands.getstatusoutput(ls)
files = output.splitlines()

# Make sure that the directory which will contain the specific datacards exists
scenarios = "scenarios"
if not os.path.exists(scenarios):
    os.popen('mkdir "%s"' %scenarios)

# Make sure that the directory which will contain the limit scan results exists
results = "results"
if not os.path.exists(results):
    os.popen('mkdir "%s"' %results)

# Do the scan!
for file in files:
    scan(file, blind)

commands.getstatusoutput("rm roostats-*")
commands.getstatusoutput("rm higgsCombineTest.Asymptotic.mH120.root")
