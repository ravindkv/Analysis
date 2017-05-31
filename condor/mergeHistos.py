
#//////////////////////////////////////////////////
#                                                 #
# Copy Ntuple root files from T2_IN_TIFR to T3.   #
# Merge them into a single root file at T3.       #
# Send the merged file back to T2_IN_TIFR.        #
# Store the full path of merged file for Analysis #
#                                                 #
#//////////////////////////////////////////////////

import os
import sys
import datetime

#USERS INPUTS
isData = True
isMC = True

mc = [
    "TTJets",
    "HplusM120",
    "ST_tW",
    "ST_t_",
    "ST_s",
    "WJetsToLNu",
    "W1JetsToLNu",
    "W2JetsToLNu",
    "W3JetsToLNu",
    "W4JetsToLNu",
    "DYJetsToLL",
    "DY1JetsToLL",
    "DY2JetsToLL",
    "DY3JetsToLL",
    "DY4JetsToLL",
    "WW",
    "WZ",
    "ZZ",
    "QCD_Pt-15to20_Mu",
    "QCD_Pt-20to30_Mu",
    "QCD_Pt-30to50_Mu",
    "QCD_Pt-50to80_Mu",
    "QCD_Pt-80to120_Mu",
    "QCD_Pt-120to170_Mu",
    "QCD_Pt-170to300_Mu",
    "HplusM80",
    "HplusM90",
    "HplusM100",
    "HplusM140",
    "HplusM150",
    "HplusM155",
    "HplusM160"]

data = [
       "MuRunB2v2",
       "MuRunCv1",
       "MuRunDv1",
       "MuRunEv1",
       "MuRunFv1",
       "MuRunGv1",
       "MuRunH2v1",
       "MuRunH3v1"
        ]

#-------------------------------
def execme(command):
    print ""
    print "\033[01;32m"+ "Excecuting: "+ "\033[00m",  command
    print ""
    os.system(command)

if isData:
    for samp in range(len(data)):
        execme("hadd -k "+str(data[samp])+"_Merged.root "+str(data[samp])+"*")

if isMC:
    for samp in range(len(mc)):
        execme("hadd -k "+str(mc[samp])+"_Merged.root "+str(mc[samp])+"*")
