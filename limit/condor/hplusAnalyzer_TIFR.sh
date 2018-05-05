#!/bin/bash

#REFERENCE
#https://github.com/florez/CONDOR

#------------------------------------------------
#pass the arguments to Analysis/runMe.sh script 
#these arguments will go to the hplusAnalyzer.C
#------------------------------------------------

massPoint=$1
channel=$2
category=$3
date

#------------------------------------------------
#this script runs on some remote condor machine.
#link lxplus to this remote machine, using scram.
#copy the compiled lxplus package to this machine

#//////////// T3 /////////////////////////////////
echo "CONDOR DIR: $_CONDOR_SCRATCH_DIR"
cd ${_CONDOR_SCRATCH_DIR}
cp -r /home/rverma/t3store2/combine/CMSSW_8_0_25/ .

#------------------------------------------------
#copy the lxplus package to the remote machine
#and run the codes at remote machine
#------------------------------------------------
cd CMSSW_8_0_25/src/HiggsAnalysis/LimitBinned/
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
./runMe.sh $massPoint $channel $category

#---------------------------------------------
#copy the output from remote machine to the lxplus
#or to any other place e.g. Tier-2
#Remove the package, after copying the output
#------------------------------------------------
echo "OUTPUT: "
ls ${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_25/src/HiggsAnalysis/LimitBinned/$channel/$category
cp -rf ${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_25/src/HiggsAnalysis/LimitBinned/$channel/$category/* /home/rverma/t3store2/condor_out/limit_condor/limit_out/$channel/$category/ 

#xrdcp -f -R ${_CONDOR_SCRATCH_DIR}/CMSSW_8_0_25/src/Analysis/13TeV/$category root://se01.indiacms.res.in:1094//cms/store/user/rverma/histo_MuMC_MuData_20170608_TIFR/
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_8_0_25

echo "DONE"
date
