#!/bin/bash
#REFERENCE
#https://github.com/florez/CONDOR

#//////////////////////////////////////
#                                     #
# hplusCond.cmd runs hplusAnalyzer.sh #
# hplusAnalyzer.sh runs runMe.sh      #
# runMe.sh runs hplusAnalyzer.C       #
#                                     #
#//////////////////////////////////////

#------------------------------------------------
#create a directory where all the outputs will be
#stored, for different merged ntuple input files
#------------------------------------------------

mkdir "outcond"
cp hplusCond.cmd "outcond"
cp hplusAnalyzer.sh "outcond"
cp mergedNtupleT2.txt "outcond"
cd "outcond"

#------------------------------------------------
#read the file, where paths of ntuples are stored
#do not put empty lines in mergedNtupleT2.txt
#------------------------------------------------

cat mergedNtupleT2.txt | while read ntupleT2Path
do
  #----------------------------------------------
  #print T2Paths of ntuple, on terminal
  #split the T2Paths of ntuples, into an array
  #get the  second last entry of the array
  #remove .root, from the input ntuple
  #----------------------------------------------
  
  echo -e "\033[01;32m input ntuple: \033[00m" $ntupleT2Path
  IFS='/ ' read -r -a array <<< "$ntupleT2Path"
  len=${#array[@]}
  sec_last=`expr $len - 1`
  ntuple=${array[$sec_last]}
  iFile=${ntuple/.root/""}
 
  #----------------------------------------------
  #copy condor scripts to each input ntuple dir
  #replace hplusCond.cmd arguments, as per input
  #submit the condor jobs, for each ntuple
  #----------------------------------------------
  
  mkdir $iFile
  cp hplusCond.cmd $iFile
  cp hplusAnalyzer.sh $iFile
  cd $iFile 
  sed -i "s:FNAME:$ntupleT2Path:g" hplusCond.cmd
  sed -i "s:OUTPUTFILE:$iFile:g" hplusCond.cmd
  sed -i "s:OUTPUTDIR:$iFile:g" hplusCond.cmd
  condor_submit hplusCond.cmd
  cd ../
done
