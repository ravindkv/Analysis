inNtupleFile=$1
outAnalFile=$2
outAnalDir=$3
mkdir $outAnalDir

#sed -i "s:outFile_.root:$fName:g" hplusAnalyzer.C
sed -i "s:inNtupleFile:$inNtupleFile:g" hplusAnalyzer.C
sed -i "s:outAnalFile:$outAnalFile:g" hplusAnalyzer.C
sed -i "s:outAnalDir:$outAnalDir:g" hplusAnalyzer.C
root -l 'runMe.C("hplusAnalyzer")'
