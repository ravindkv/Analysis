#!/bin/bash
#REFERENCE
#https://github.com/florez/CONDOR

#//////////////////////////////////////
#                                     #
# hplusCond.sub runs hplusAnalyzer_TIFR.sh #
# hplusAnalyzer_TIFR.sh runs runMe.sh      #
# runMe.sh runs hplusAnalyzer.C       #
#                                     #
#//////////////////////////////////////

#------------------------------------------------
#create a directory where all the outputs will be
#stored, for different merged ntuple input files
#------------------------------------------------
ntupleT2Paths=$1

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

log="CondorLog_"
logDir=$log${ntupleT2Paths/.txt/""}
baseDir="/home/rverma/t3store/AN-18-061/CondorLog/MLFitCondorLog/"
mkdir -p $baseDir$logDir
outcond="$baseDir$logDir"

cp hplusCond.sub $outcond
cp hplusAnalyzer_TIFR.sh $outcond
cp $ntupleT2Paths $outcond
cd $outcond

#------------------------------------------------
#read the file, where paths of ntuples are stored
#do not put empty lines in mergedNtupleT2.txt
#------------------------------------------------
count=0
cat $ntupleT2Paths | while read ntupleT2Path
do
  #----------------------------------------------
  #print T2Paths of ntuple, on terminal
  #split the T2Paths of ntuples, into an array
  #get the  second last entry of the array
  #remove .root, from the input ntuple
  #----------------------------------------------
  ((count++))
  echo -e "\033[01;32m input ntuple=\033[00m" $count": " $ntupleT2Path
  IFS='/' read -r -a array <<< "$ntupleT2Path"
  len=${#array[@]}
  limitDir="${array[0]}"
  channel="${array[1]}"
  cat="${array[2]}"
  mass="${array[3]}"

  mkdir -p /home/rverma/t3store/AN-18-061/CondorOut/MLFitCondorOut/$limitDir/$channel/$cat
  #echo $ntuple
  iFile=$limitDir$channel$cat$mass
 
  #----------------------------------------------
  #copy condor scripts to each input ntuple dir
  #replace hplusCond.sub arguments, as per input
  #submit the condor jobs, for each ntuple
  #----------------------------------------------
  mkdir -p $iFile
  cp hplusCond.sub $iFile
  cp hplusAnalyzer_TIFR.sh $iFile
  cd $iFile 
  sed -i "s:LIMITDIR:$limitDir:g" hplusCond.sub
  sed -i "s:CHANNEL:$channel:g" hplusCond.sub
  sed -i "s:CAT:$cat:g" hplusCond.sub
  sed -i "s:MASS:$mass:g" hplusCond.sub
  condor_submit hplusCond.sub
  cd ../
done
