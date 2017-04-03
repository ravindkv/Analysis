#!/bin/bash
#for iFile in 'mergedNtupleT2.txt'; 
mkdir "outcond"
cp hplusCond.cmd "outcond"
cp hplusAnalyzer.sh "outcond"
cp mergedNtupleT2.txt "outcond"
cd "outcond"

cat mergedNtupleT2.txt | while read ntuple
do
  echo $ntuple
  iFile=${ntuple/.root/""}
  
  mkdir $iFile
  cp hplusCond.cmd hplusCond_"$iFile".cmd
  mv hplusCond_"$iFile".cmd $iFile
  cp hplusAnalyzer.sh $iFile
  
  cd $iFile 
  outFile="analHistos_$iFile"
  sed -i "s:FNAME:$ntuple:g" hplusCond_"$iFile".cmd
  sed -i "s:OUTPUTFILE:$outFile:g" hplusCond_"$iFile".cmd
  sed -i "s:OUTPUTDIR:$iFile:g" hplusCond_"$iFile".cmd
  condor_submit hplusCond_"$iFile".cmd
  cd ../
done

