#!/bin/bash
#https://github.com/florez/CONDOR

# These input parameters are going to be passed to the 
# script from the submit_condor_jobs.cmd script, 
# which takes the input from the run_code.sh 
 
inNtupleFile=$1
outAnalFile=$2
outAnalDir=$3
date

#cd to the are where you have installed your CMSSW release, e.g:
cd /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src
#cd ../../../src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

#Copy the directory where you have the compiled code, e.g.:
cp -r /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src/Analysis .
#cp -r Analysis .
cd Analysis
./runMe.sh $inNtupleFile $outAnalFile $outAnalDir

