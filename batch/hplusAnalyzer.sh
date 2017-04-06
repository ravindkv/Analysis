#!/bin/bash
#https://github.com/florez/CONDOR

inNtupleFile=$1
outAnalFile=$2
outAnalDir=$3
date

#ls $PWD
#cd to the are where you have installed your CMSSW release, e.g:
cd /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src
#cmsenv
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
echo "CONDOR DIR: $_CONDOR_SCRATCH_DIR"
cd ${_CONDOR_SCRATCH_DIR}

echo "CONDOR SCRATCH DIR: $PWD"
cp -r /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src/Analysis .

#Copy the directory where you have the compiled code, e.g.:
#cp -r /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src/Analysis .
cd Analysis
./runMe.sh $inNtupleFile $outAnalFile $outAnalDir

echo "LIST BEFORE MOVING"
ls ${_CONDOR_SCRATCH_DIR}/Analysis
ls ${_CONDOR_SCRATCH_DIR}/Analysis/13TeV

xrdcp  ${_CONDOR_SCRATCH_DIR}/Analysis/13TeV/* /afs/cern.ch/work/r/rverma/private/analysis/CMSSW_7_2_3/src
