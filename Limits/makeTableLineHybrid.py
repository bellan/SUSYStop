#! /usr/bin/env python


import sys, os, commands, math, ROOT, tokenize


def getLimitValue(filename):
    file = open(filename,'r')
    for line in file:
        # find the line with the input

        found =  line.find("@ 95% CL")
        if found >= 0:
            # split the line using comma as separator
            spline =  line.split(" ")
            return float(spline[3])
            

def printLine(name,array):
    print "{0:s}:     & {1:.1f} & ${2:.1f}_".format(name,array["Obs"],array["Exp"])+"{"+"{0:.1f}".format(array["Expm1s"]-array["Exp"])+"}^{"+"+ {0:.1f}".format(array["Expp1s"]-array["Exp"])+"}$ \\"+"{0:s}".format("\\")
    


def printRegion(name,filebasename):

    array = {'Obs': -1,'Exp': -1,'Expm1s': -1,'Expp1s': -1}

    array["Obs"]    = getLimitValue(filebasename+"_hybr_obs.txt")
    array["Exp"]    = getLimitValue(filebasename+"_hybr_exp.txt")
    array["Expm1s"] = getLimitValue(filebasename+"_hybr_exp_m1s.txt")
    array["Expp1s"] = getLimitValue(filebasename+"_hybr_exp_p1s.txt")

    printLine(name, array)




printRegion("Baseline","baseline")
printRegion("Loose","loose")
printRegion("Medium","medium")
printRegion("Tight","tight")
