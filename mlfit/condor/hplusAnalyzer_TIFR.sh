#!/bin/bash

#REFERENCE
#https://github.com/florez/CONDOR

#------------------------------------------------
#pass the arguments to Analysis/runMe.sh script 
#these arguments will go to the hplusAnalyzer.C
#------------------------------------------------

limitDir=$1
channel=$2
category=$3
massPoint=$4
date

#------------------------------------------------
#this script runs on some remote condor machine.
#link lxplus to this remote machine, using scram.
#copy the compiled lxplus package to this machine

#//////////// T3 /////////////////////////////////
echo "CONDOR DIR: $_CONDOR_SCRATCH_DIR"
cd ${_CONDOR_SCRATCH_DIR}
cp -r /home/rverma/t3store/AN-18-061/ExclusionLimit/CMSSW_8_0_26/ .

#------------------------------------------------
#copy the lxplus package to the remote machine
#and run the codes at remote machine
#------------------------------------------------
cd CMSSW_8_0_26/src/HiggsAnalysis/HplusTocs13TeVMLFit/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
python MyLiklihoodFit.py --limit $limitDir --channel $channel --hist $category --mass $massPoint

#---------------------------------------------
#copy the output from remote machine to the lxplus
#or to any other place e.g. Tier-2
#Remove the package, after copying the output
#------------------------------------------------
echo "OUTPUT: "
ls ${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_26/src/HiggsAnalysis/HplusTocs13TeVMLFit/$limitDir/$channel/$category
cp -rf ${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_26/src/HiggsAnalysis/HplusTocs13TeVMLFit/$limitDir/$channel/$category/* /home/rverma/t3store/AN-18-061/CondorOut/MLFitCondorOut/$limitDir/$channel/$category/ 

cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_8_0_26

echo "DONE"
date

