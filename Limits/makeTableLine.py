#! /usr/bin/env python


import sys, os, commands, math, ROOT, tokenize



def printLine(filename,name):
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


    print "{0:s}:     & {1:.1f} & ${2:.1f}_".format(name,array["Obs"],array["Exp"])+"{"+"{0:.1f}".format(array["Expm1s"]-array["Exp"])+"}^{"+"+ {0:.1f}".format(array["Expp1s"]-array["Exp"])+"}$ \\"+"{0:s}".format("\\")


list = [0,5,10,15,20,25,30]

for i in list:
    print "\underline{Asymptotic CLs ("+"{0:.1f}".format(i)+" \%):} & & \\"+"{0:s}".format("\\")
    printLine("results/baseline_"+str(i)+".txt","Baseline")
    #printLine("results/loose_"+str(i)+".txt","Loose")
    #printLine("results/medium_"+str(i)+".txt","Medium")
    #printLine("results/right_"+str(i)+".txt","Tight")
