#!/bin/bash
#http://research.cs.wisc.edu/htcondor/manual/index.html
mkdir "outcond"
cp hplusCond.cmd "outcond"
cp hplusAnalyzer.sh "outcond"
cp mergedNtupleT2.txt "outcond"
cd "outcond"
cat mergedNtupleT2.txt | while read ntuple
do
  echo "\033[01;32m input Ntuple: \033[00m" $ntuple
  iFile=${ntuple/.root/""}
  
  mkdir $iFile
  cp hplusCond.cmd $iFile
  cp hplusAnalyzer.sh $iFile
  
  cd $iFile 
  outFile="analHistos_$iFile"
  sed -i "s:FNAME:$ntuple:g" hplusCond.cmd
  sed -i "s:OUTPUTFILE:$outFile:g" hplusCond.cmd
  sed -i "s:OUTPUTDIR:$iFile:g" hplusCond.cmd
  condor_submit hplusCond.cmd
  cd ../
done

