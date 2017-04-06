inNtupleFile=$1
outAnalFile=$2
outAnalDir=$3
cd 13TeV
mkdir $outAnalDir
cd ../

echo "INPUT FILE: " $inNtupleFile
cp hplusAnalyzer.C ../
sed -i "s:inputFile:$inNtupleFile:g" hplusAnalyzer.C
sed -i "s:outputFile:$outAnalFile:g" hplusAnalyzer.C
sed -i "s:outputDir:$outAnalDir:g" hplusAnalyzer.C
root -l 'runMe.C("hplusAnalyzer")'
echo "SUCCESSFULL"
mv ../hplusAnalyzer.C .
